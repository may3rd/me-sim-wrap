# Phase 2B Process Units Design

## Decision

Extend `src/mesim/units.py`. Do not import the old `unit-converter`, port DWSIM conversion code, or add an arbitrary unit-expression parser.

## Boundary

The registry accepts exact approved symbols and aliases, validates semantic dimensions, converts to one documented SI base, and preserves the submitted value and unit in `Quantity`. Unknown expressions fail instead of being guessed.

Affine units remain explicit:

- absolute temperature uses offsets
- temperature difference never uses offsets
- gauge-pressure symbols use the fixed standard-atmosphere reference and convert to absolute Pa

Phase 2B covers process dimensions needed through the internal-alpha roadmap. Later equipment adds units only with a failing calculation-boundary test.

## Reuse

Reuse only independently verified factors from `packages/unit-converter/unit_converter/data.py`. Its parser and unit arithmetic are rejected because division and affine compound units are incorrect.

## Verification

Every added dimension gets one SI and one non-SI equation vector, a round trip, a wrong-dimension rejection, and non-finite rejection through the shared boundary.
