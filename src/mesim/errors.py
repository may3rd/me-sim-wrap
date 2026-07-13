class ValidationError(ValueError):
    """Input or model data violates a calculation contract."""


class MissingCompoundData(ValidationError):
    """A required compound property is unavailable."""


class OutOfRangeError(ValidationError):
    """A model input is outside its declared validity range."""


class ConvergenceError(RuntimeError):
    """An iterative calculation did not converge."""
