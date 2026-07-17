from collections.abc import Callable, Mapping
from dataclasses import dataclass
import math
from types import MappingProxyType

from .errors import ConvergenceError, ValidationError


PortValues = Mapping[str, object]


def _ports(value: tuple[str, ...], name: str) -> None:
    if len(set(value)) != len(value) or any(not isinstance(port, str) or not port for port in value):
        raise ValidationError(f"{name} must be unique non-empty strings")


@dataclass(frozen=True, slots=True)
class UnitOperation:
    tag: str
    input_ports: tuple[str, ...]
    output_ports: tuple[str, ...]
    calculate: Callable[[PortValues], PortValues]

    def __post_init__(self) -> None:
        if not isinstance(self.tag, str) or not self.tag:
            raise ValidationError("unit tag must be a non-empty string")
        _ports(self.input_ports, "unit input ports")
        _ports(self.output_ports, "unit output ports")
        if not self.output_ports or not callable(self.calculate):
            raise ValidationError("unit requires output ports and a calculation")


@dataclass(frozen=True, slots=True)
class Connection:
    source_tag: str
    source_port: str
    target_tag: str
    target_port: str


@dataclass(frozen=True, slots=True)
class Flowsheet:
    units: tuple[UnitOperation, ...]
    connections: tuple[Connection, ...]


@dataclass(frozen=True, slots=True)
class FlowsheetResult:
    values: Mapping[tuple[str, str], object]


@dataclass(frozen=True, slots=True)
class RecycleIteration:
    iteration: int
    guess: tuple[float, ...]
    calculated: tuple[float, ...]
    residual: tuple[float, ...]
    scaled_norm: float
    damping: float


@dataclass(frozen=True, slots=True)
class RecycleResult:
    values: tuple[float, ...]
    history: tuple[RecycleIteration, ...]
    algorithm: str = "direct_substitution"


class RecycleConvergenceError(ConvergenceError):
    def __init__(self, message: str, history: tuple[RecycleIteration, ...]):
        super().__init__(message)
        self.history = history


def solve_recycle(
    evaluate: Callable[[tuple[float, ...]], tuple[float, ...]],
    initial_guess: tuple[float, ...],
    scales: tuple[float, ...],
    tolerances: tuple[float, ...],
    damping: float = 1.0,
    max_iterations: int = 100,
) -> RecycleResult:
    """Solve explicit tear variables by damped direct substitution with full history."""
    if not callable(evaluate):
        raise ValidationError("recycle evaluation must be callable")
    try:
        guess = tuple(initial_guess)
        scale_values = tuple(scales)
        tolerance_values = tuple(tolerances)
    except TypeError as error:
        raise ValidationError("recycle vectors must be finite sequences") from error
    if not guess or len(guess) != len(scale_values) or len(guess) != len(tolerance_values):
        raise ValidationError("recycle guess, scales, and tolerances must have the same non-zero length")
    if any(
        isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
        for value in guess
    ) or any(
        isinstance(value, bool) or not isinstance(value, (int, float))
        or not math.isfinite(value) or value <= 0.0
        for value in scale_values + tolerance_values
    ):
        raise ValidationError("recycle vectors must contain finite values and positive scales and tolerances")
    if (
        isinstance(damping, bool) or not isinstance(damping, (int, float))
        or not math.isfinite(damping) or not 0.0 < damping <= 1.0
        or isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations <= 0
    ):
        raise ValidationError("recycle damping and maximum iterations are invalid")

    guess = tuple(float(value) for value in guess)
    scale_values = tuple(float(value) for value in scale_values)
    tolerance_values = tuple(float(value) for value in tolerance_values)
    history: list[RecycleIteration] = []
    for iteration in range(1, max_iterations + 1):
        try:
            calculated = tuple(evaluate(guess))
        except TypeError as error:
            raise ValidationError("recycle evaluation must return a finite sequence") from error
        if len(calculated) != len(guess) or any(
            isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
            for value in calculated
        ):
            raise ValidationError("recycle evaluation returned an invalid vector")
        calculated = tuple(float(value) for value in calculated)
        residual = tuple(value - prior for value, prior in zip(calculated, guess))
        scaled_norm = max(abs(value) / scale for value, scale in zip(residual, scale_values))
        history.append(RecycleIteration(
            iteration, guess, calculated, residual, scaled_norm, float(damping),
        ))
        if all(abs(value) <= tolerance for value, tolerance in zip(residual, tolerance_values)):
            return RecycleResult(calculated, tuple(history))
        guess = tuple(prior + damping * value for prior, value in zip(guess, residual))

    raise RecycleConvergenceError(
        "recycle direct substitution did not converge", tuple(history),
    )


def solve_energy_recycle(
    evaluate_duty_w: Callable[[float], float],
    initial_duty_w: float,
    scale_w: float,
    tolerance_w: float,
    damping: float = 1.0,
    max_iterations: int = 100,
) -> RecycleResult:
    """Solve a scalar energy-stream tear in watts by direct substitution."""
    if not callable(evaluate_duty_w):
        raise ValidationError("energy recycle evaluation must be callable")

    def evaluate(values: tuple[float, ...]) -> tuple[float, ...]:
        return (evaluate_duty_w(values[0]),)

    return solve_recycle(
        evaluate,
        initial_guess=(initial_duty_w,),
        scales=(scale_w,),
        tolerances=(tolerance_w,),
        damping=damping,
        max_iterations=max_iterations,
    )


def _validate(flowsheet: Flowsheet) -> dict[str, UnitOperation]:
    units = {unit.tag: unit for unit in flowsheet.units}
    if not units or len(units) != len(flowsheet.units):
        raise ValidationError("flowsheet unit tags must be unique and non-empty")
    targets: set[tuple[str, str]] = set()
    sources: set[tuple[str, str]] = set()
    for connection in flowsheet.connections:
        try:
            source, target = units[connection.source_tag], units[connection.target_tag]
        except KeyError as error:
            raise ValidationError("connection references an unknown unit") from error
        source_key = (connection.source_tag, connection.source_port)
        target_key = (connection.target_tag, connection.target_port)
        if connection.source_port not in source.output_ports or connection.target_port not in target.input_ports:
            raise ValidationError("connection references an unknown port")
        if source_key in sources or target_key in targets:
            raise ValidationError("each flowsheet port may have only one connection")
        sources.add(source_key)
        targets.add(target_key)
    if any((unit.tag, port) not in targets for unit in units.values() for port in unit.input_ports):
        raise ValidationError("every unit input port requires one connection")
    return units


def topological_order(flowsheet: Flowsheet) -> tuple[str, ...]:
    units = _validate(flowsheet)
    dependencies = {tag: set() for tag in units}
    successors = {tag: set() for tag in units}
    for connection in flowsheet.connections:
        dependencies[connection.target_tag].add(connection.source_tag)
        successors[connection.source_tag].add(connection.target_tag)
    ready = sorted(tag for tag, dependency in dependencies.items() if not dependency)
    order: list[str] = []
    while ready:
        tag = ready.pop(0)
        order.append(tag)
        for successor in sorted(successors[tag]):
            dependencies[successor].remove(tag)
            if not dependencies[successor]:
                ready.append(successor)
        ready.sort()
    if len(order) != len(units):
        raise ValidationError("flowsheet cycles are unsupported")
    return tuple(order)


def solve(flowsheet: Flowsheet) -> FlowsheetResult:
    units = _validate(flowsheet)
    order = topological_order(flowsheet)
    inputs = {(connection.target_tag, connection.target_port): (connection.source_tag, connection.source_port) for connection in flowsheet.connections}
    values: dict[tuple[str, str], object] = {}
    for tag in order:
        unit = units[tag]
        supplied = MappingProxyType({port: values[inputs[(tag, port)]] for port in unit.input_ports})
        calculated = unit.calculate(supplied)
        if not isinstance(calculated, Mapping) or set(calculated) != set(unit.output_ports):
            raise ValidationError(f"unit {tag} must return exactly its declared output ports")
        values.update({(tag, port): calculated[port] for port in unit.output_ports})
    return FlowsheetResult(MappingProxyType(values))
