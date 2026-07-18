"""DWSIM/ChemSep temperature-equation evaluation shared by property families."""

from __future__ import annotations

import math

from ..errors import ValidationError


SUPPORTED_EQUATIONS = frozenset((1, 2, 3, 4, 10, 16, 100, 101, 102, 105, 106, 116))


def evaluate_temperature_equation(
    equation: int,
    coefficients: tuple[float, float, float, float, float],
    temperature_k: float,
    critical_temperature_k: float | None = None,
) -> float:
    """Evaluate the matching DWSIM ``CalcCSTDepProp`` equation."""
    if equation not in SUPPORTED_EQUATIONS:
        raise ValidationError(f"unsupported temperature equation {equation}")
    if not math.isfinite(temperature_k) or temperature_k <= 0:
        raise ValidationError("absolute temperature must be finite and positive")
    a, b, c, d, e = coefficients
    if equation in (105, 106, 116):
        if critical_temperature_k is None or not math.isfinite(critical_temperature_k) or critical_temperature_k <= 0:
            raise ValidationError("a positive finite critical temperature is required")
    reduced = None if critical_temperature_k is None else temperature_k / critical_temperature_k
    if equation == 105 and (b <= 0 or c <= 0 or temperature_k >= c):
        raise ValidationError("equation 105 state is outside its real-valued domain")
    if equation in (106, 116) and reduced >= 1.0:
        raise ValidationError(
            f"equation {equation} state is at or above the critical temperature"
        )
    try:
        if equation == 1:
            value = a
        elif equation == 2:
            value = a + b * temperature_k
        elif equation == 3:
            value = a + b * temperature_k + c * temperature_k**2
        elif equation == 4:
            value = a + b * temperature_k + c * temperature_k**2 + d * temperature_k**3
        elif equation == 10:
            value = math.exp(a - b / (temperature_k + c))
        elif equation == 16:
            value = a + math.exp(
                b / temperature_k + c + d * temperature_k + e * temperature_k**2
            )
        elif equation == 100:
            value = (
                a
                + b * temperature_k
                + c * temperature_k**2
                + d * temperature_k**3
                + e * temperature_k**4
            )
        elif equation == 101:
            value = math.exp(
                a
                + b / temperature_k
                + c * math.log(temperature_k)
                + d * temperature_k**e
            )
        elif equation == 102:
            value = a * temperature_k**b / (
                1.0 + c / temperature_k + d / temperature_k**2
            )
        elif equation == 105:
            value = a / b ** (1.0 + (1.0 - temperature_k / c) ** d)
        elif equation == 106:
            value = a * (1.0 - reduced) ** (
                b + c * reduced + d * reduced**2 + e * reduced**3
            )
        else:  # 116
            one_minus = 1.0 - reduced
            value = (
                a
                + b * one_minus**0.35
                + c * one_minus ** (2.0 / 3.0)
                + d * one_minus
                + e * one_minus ** (4.0 / 3.0)
            )
    except (OverflowError, ValueError, ZeroDivisionError) as error:
        raise ValidationError("temperature equation is outside the representable range") from error
    if not math.isfinite(value):
        raise ValidationError("temperature equation result must be finite")
    return value
