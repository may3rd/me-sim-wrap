import ast
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


@dataclass(frozen=True, slots=True)
class AdjustIteration:
    iteration: int
    manipulated: float
    controlled: float
    target: float
    residual: float
    scaled_norm: float
    derivative: float | None
    step: float
    damping: float


@dataclass(frozen=True, slots=True)
class AdjustResult:
    manipulated: float
    controlled: float
    history: tuple[AdjustIteration, ...]
    algorithm: str = "bounded_newton"


class AdjustConvergenceError(ConvergenceError):
    def __init__(self, message: str, history: tuple[AdjustIteration, ...]):
        super().__init__(message)
        self.history = history


@dataclass(frozen=True, slots=True)
class SpecificationResult:
    expression: str
    source_value: float
    target_before: float
    unconstrained_value: float
    target_value: float
    minimum_value: float | None
    maximum_value: float | None
    clamped: bool


_SPEC_FUNCTIONS: Mapping[str, Callable[..., float]] = MappingProxyType({
    "abs": abs,
    "acos": math.acos,
    "asin": math.asin,
    "atan": math.atan,
    "atan2": math.atan2,
    "ceiling": math.ceil,
    "cos": math.cos,
    "cosh": math.cosh,
    "exp": math.exp,
    "floor": math.floor,
    "log": math.log,
    "log10": math.log10,
    "max": max,
    "min": min,
    "pow": math.pow,
    "round": round,
    "sin": math.sin,
    "sinh": math.sinh,
    "sqrt": math.sqrt,
    "tan": math.tan,
    "tanh": math.tanh,
})


def _evaluate_spec_node(node: ast.AST, variables: Mapping[str, float]) -> float:
    if isinstance(node, ast.Constant):
        if isinstance(node.value, bool) or not isinstance(node.value, (int, float)):
            raise ValidationError("specification constants must be numeric")
        return float(node.value)
    if isinstance(node, ast.Name):
        try:
            return variables[node.id.lower()]
        except KeyError as error:
            raise ValidationError(f"unsupported specification name: {node.id}") from error
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        value = _evaluate_spec_node(node.operand, variables)
        return value if isinstance(node.op, ast.UAdd) else -value
    if isinstance(node, ast.BinOp):
        left = _evaluate_spec_node(node.left, variables)
        right = _evaluate_spec_node(node.right, variables)
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
        if isinstance(node.op, ast.Mult):
            return left * right
        if isinstance(node.op, ast.Div):
            return left / right
        if isinstance(node.op, (ast.Pow, ast.BitXor)):
            return left**right
        if isinstance(node.op, ast.Mod):
            return left % right
        raise ValidationError("unsupported specification operator")
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and not node.keywords:
        try:
            function = _SPEC_FUNCTIONS[node.func.id.lower()]
        except KeyError as error:
            raise ValidationError(
                f"unsupported specification function: {node.func.id}"
            ) from error
        return float(function(*(_evaluate_spec_node(argument, variables) for argument in node.args)))
    raise ValidationError("unsupported specification expression")


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


def solve_adjust(
    evaluate_controlled: Callable[[float], float],
    target: float,
    initial_guess: float,
    lower_bound: float,
    upper_bound: float,
    controlled_scale: float,
    tolerance: float,
    step_size: float,
    damping: float = 1.0,
    max_iterations: int = 25,
) -> AdjustResult:
    """Solve one bounded manipulated variable with finite-difference Newton steps."""
    values = (
        target, initial_guess, lower_bound, upper_bound,
        controlled_scale, tolerance, step_size, damping,
    )
    if not callable(evaluate_controlled):
        raise ValidationError("adjust evaluation must be callable")
    if any(
        isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
        for value in values
    ) or not lower_bound < upper_bound or not lower_bound <= initial_guess <= upper_bound:
        raise ValidationError("adjust values and bounds are invalid")
    if (
        controlled_scale <= 0.0 or tolerance <= 0.0 or step_size <= 0.0
        or not 0.0 < damping <= 1.0
        or isinstance(max_iterations, bool) or not isinstance(max_iterations, int)
        or max_iterations <= 0
    ):
        raise ValidationError("adjust scales, tolerance, damping, and iterations are invalid")

    target = float(target)
    manipulated = float(initial_guess)
    controlled_scale = float(controlled_scale)
    tolerance = float(tolerance)
    step_size = float(step_size)
    damping = float(damping)
    history: list[AdjustIteration] = []

    def evaluate(value: float) -> float:
        controlled = evaluate_controlled(value)
        if (
            isinstance(controlled, bool) or not isinstance(controlled, (int, float))
            or not math.isfinite(controlled)
        ):
            raise ValidationError("adjust evaluation returned an invalid controlled value")
        return float(controlled)

    for iteration in range(1, max_iterations + 1):
        controlled = evaluate(manipulated)
        residual = target - controlled
        scaled_norm = abs(residual) / controlled_scale
        if abs(residual) <= tolerance:
            history.append(AdjustIteration(
                iteration, manipulated, controlled, target, residual,
                scaled_norm, None, 0.0, damping,
            ))
            return AdjustResult(manipulated, controlled, tuple(history))

        low_probe = max(float(lower_bound), manipulated - step_size)
        high_probe = min(float(upper_bound), manipulated + step_size)
        if high_probe <= low_probe:
            history.append(AdjustIteration(
                iteration, manipulated, controlled, target, residual,
                scaled_norm, None, 0.0, damping,
            ))
            raise AdjustConvergenceError("adjust has no finite-difference interval", tuple(history))
        derivative = (evaluate(high_probe) - evaluate(low_probe)) / (high_probe - low_probe)
        if not math.isfinite(derivative) or derivative == 0.0:
            history.append(AdjustIteration(
                iteration, manipulated, controlled, target, residual,
                scaled_norm, derivative, 0.0, damping,
            ))
            raise AdjustConvergenceError("adjust derivative is singular", tuple(history))

        proposed = manipulated + damping * residual / derivative
        bounded = min(float(upper_bound), max(float(lower_bound), proposed))
        applied_step = bounded - manipulated
        history.append(AdjustIteration(
            iteration, manipulated, controlled, target, residual,
            scaled_norm, derivative, applied_step, damping,
        ))
        if applied_step == 0.0:
            raise AdjustConvergenceError("adjust is pinned at a bound", tuple(history))
        manipulated = bounded

    raise AdjustConvergenceError("adjust bounded Newton solve did not converge", tuple(history))


def apply_specification(
    expression: str,
    source_value: float,
    target_value: float,
    minimum_value: float | None = None,
    maximum_value: float | None = None,
) -> SpecificationResult:
    """Evaluate a DWSIM-style X/Y scalar expression and apply optional bounds."""
    if not isinstance(expression, str) or not expression.strip():
        raise ValidationError("specification expression must be a non-empty string")
    if len(expression) > 1000:
        raise ValidationError("specification expression is too long")
    if any(
        isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
        for value in (source_value, target_value)
    ):
        raise ValidationError("specification source and target values must be finite numbers")
    if (minimum_value is None) != (maximum_value is None):
        raise ValidationError("specification bounds must be supplied together")
    if minimum_value is not None and (
        isinstance(minimum_value, bool) or not isinstance(minimum_value, (int, float))
        or isinstance(maximum_value, bool) or not isinstance(maximum_value, (int, float))
        or not math.isfinite(minimum_value) or not math.isfinite(maximum_value)
        or minimum_value > maximum_value
    ):
        raise ValidationError("specification bounds are invalid")

    try:
        parsed = ast.parse(expression, mode="eval")
    except (SyntaxError, ValueError) as error:
        raise ValidationError("specification expression is invalid") from error
    if sum(1 for _ in ast.walk(parsed)) > 64:
        raise ValidationError("specification expression is too complex")
    variables = MappingProxyType({
        "x": float(source_value),
        "y": float(target_value),
        "pi": math.pi,
        "e": math.e,
    })
    try:
        unconstrained = _evaluate_spec_node(parsed.body, variables)
    except ValidationError:
        raise
    except (ArithmeticError, OverflowError, TypeError, ValueError) as error:
        raise ValidationError("specification expression could not be evaluated") from error
    if not math.isfinite(unconstrained):
        raise ValidationError("specification expression returned a non-finite value")

    minimum = None if minimum_value is None else float(minimum_value)
    maximum = None if maximum_value is None else float(maximum_value)
    calculated = unconstrained
    if minimum is not None and maximum is not None:
        calculated = min(maximum, max(minimum, unconstrained))
    return SpecificationResult(
        expression=expression,
        source_value=float(source_value),
        target_before=float(target_value),
        unconstrained_value=unconstrained,
        target_value=calculated,
        minimum_value=minimum,
        maximum_value=maximum,
        clamped=calculated != unconstrained,
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
