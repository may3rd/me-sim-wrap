"""PR/Lee-Kesler equilibrium boundary used by DWSIM's hybrid package."""

from ..compounds import Compound, PRInteractions
from .flash import TPFlashResult, tp_flash
from .peng_robinson import PRMixtureState, PengRobinsonMixture


class PengRobinsonLeeKeslerMixture(PengRobinsonMixture):
    """PR fugacity model from PR/LK; Lee-Kesler calorics are out of scope here.

    DWSIM's hybrid package inherits its phase-equilibrium calculation from the
    classic Peng-Robinson package and overrides caloric and compressibility
    property paths with Lee-Kesler corresponding states.  Keeping this type
    distinct prevents callers from assuming that classic-PR calorics have been
    accepted for the hybrid package.
    """


def pr_lk_tp_flash(
    compounds: tuple[Compound, ...],
    composition: tuple[float, ...],
    interactions: PRInteractions,
    temperature_k: float,
    pressure_pa: float,
    *,
    max_iterations: int = 100,
) -> TPFlashResult:
    """Run the PR equilibrium flash used by the DWSIM PR/LK package."""
    return tp_flash(
        compounds,
        composition,
        interactions,
        temperature_k,
        pressure_pa,
        max_iterations=max_iterations,
    )


__all__ = (
    "PRMixtureState",
    "PengRobinsonLeeKeslerMixture",
    "TPFlashResult",
    "pr_lk_tp_flash",
)
