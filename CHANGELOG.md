# Changelog

All notable changes are documented here. This project follows semantic versioning for the Python package, while API, data, golden-case, and solver schemas retain independent identifiers.

## [Unreleased]

- Await hosted Linux amd64, Linux arm64, and macOS arm64 release-workflow results before promoting the candidate or creating the `0.1.0` tag.

## [0.1.0rc1] - 2026-07-18

### Added

- Explicit SI unit boundary, versioned compound/correlation/interaction/reaction data, Peng-Robinson thermodynamics, TP/PH flashes, and immutable streams.
- Golden-backed basic operations, pressure changers, heat exchangers, hydraulics, reactors, recycles/logical blocks, columns, dynamics/PID, and solar/wind/hydroelectric source-equation slices.
- FastAPI U0 boundary, request-size and calculation-deadline guards, deterministic Windows DWSIM capture tooling, cross-platform release workflow, and representative performance benchmark.
- Release support matrix, model limitations, and data-provenance records.

### Changed

- Replaced per-request interpreter startup with a reusable two-worker spawned process pool while preserving hard timeout termination.

### Known differences

- Bubble/dew pressure lacks a valid nonzero DWSIM property capture.
- Several capabilities are partial fixed-thermodynamic or source-equation gates rather than general predictive models; see `docs/model-limitations.md`.
- Public distribution remains subject to GPLv3-derived-work review.
