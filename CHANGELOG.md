# Changelog

All notable changes are documented here. This project follows semantic versioning for the Python package, while API, data, golden-case, and solver schemas retain independent identifiers.

## [Unreleased]

### Added

- Deterministic, nonzero DWSIM bubble/dew-pressure captures and executable methane/ethane PR parity gates.
- A saved-source acetone/methanol NRTL activity and modified-Raoult bubble/dew slice for the accepted pressure-swing column domain.
- Live NRTL K-value and bubble-point-temperature parity across all 20 stages of the accepted acetone column profile.

## [0.1.0] - 2026-07-18

### Added

- Explicit SI unit boundary, versioned compound/correlation/interaction/reaction data, Peng-Robinson thermodynamics, TP/PH flashes, and immutable streams.
- Golden-backed basic operations, pressure changers, heat exchangers, hydraulics, reactors, recycles/logical blocks, columns, dynamics/PID, and solar/wind/hydroelectric source-equation slices.
- Methanol-carbonylation CSTR parity with deterministic, property-error-free DWSIM captures and explicit mol/kmol kinetic-unit conversion.
- FastAPI U0 boundary, request-size and calculation-deadline guards, deterministic Windows DWSIM capture tooling, cross-platform release workflow, and representative performance benchmark.
- Release support matrix, model limitations, and data-provenance records.

### Changed

- Replaced per-request interpreter startup with a reusable two-worker spawned process pool while preserving hard timeout termination.

### Known differences

- Bubble/dew pressure lacks a valid nonzero DWSIM property capture.
- Several capabilities are partial fixed-thermodynamic or source-equation gates rather than general predictive models; see `docs/model-limitations.md`.
- Linux amd64, Linux arm64, macOS arm64, and container jobs are defined but remain unexecuted because GitHub rejected the jobs for an account-level billing lock; this source promotion is not a multi-architecture verification or release tag.
- Public distribution remains subject to GPLv3-derived-work review.
