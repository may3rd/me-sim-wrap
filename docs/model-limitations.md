# Model limitations

This document is a release boundary for `me-sim-wrap` 0.1.0. A result is supported only when its model, phase domain, compound data, and solver mode are explicitly covered in [compatibility.md](compatibility.md). Unsupported input must be rejected rather than silently approximated.

## Safety and design use

- The kernel is an internal calculation component, not a certified process-design, relief-design, hazard-analysis, or operator-training product.
- API RP 520 utility equations reproduce the DWSIM source path for frozen test vectors only. They do not size a complete relief system, validate standards editions, or replace engineering review.
- Specialty renewable-energy equations are source-equation comparisons, not site-yield forecasts or manufacturer guarantees.

## Thermodynamics and compounds

- Classic Peng-Robinson is the principal full flash-and-caloric equation of state. PR78, advanced PR78, PRSV2-M, and PRSV2-VL currently expose phase states; SRK and advanced SRK expose phase states and TP flash; PR/Lee-Kesler exposes its PR equilibrium side but not Lee-Kesler calorics. The advanced cubic Mercury record is EOS-only and does not add general pure-property capabilities. Electrolytes, solids, petroleum characterization, hydrates, and CAPE-OPEN packages are not Python property packages. Activity-coefficient coverage consists of the documented acetone/methanol NRTL pressure/caloric slice, acetone/methanol Wilson activities, and 1-propanol/water UNIQUAC, UNIFAC, UNIFAC-LL, Dortmund, and NIST Modified UNIFAC activities; none is a general activity-property package or LLE split solver.
- The thermodynamic-system registry is fixed to the model IDs documented in `thermodynamic-systems.md`; it is a capability boundary, not a runtime plugin or automatic model-selection system.
- The initial five-compound PR domain is methane, ethane, propane, n-butane, and n-pentane. Additional catalog records exist to support specific golden-backed reaction and column cases; their presence does not imply universal model coverage.
- Correlation temperature bounds are authoritative. Extrapolation is rejected unless an API explicitly exposes and documents an opt-in flag.
- Unmodified cubic density and Peneloux-translated density are separate model choices. DWSIM stream density fields may use Rackett, experimental data, or volume translation and are not interchangeable.
- PR bubble/dew solvers have deterministic, nonzero DWSIM parity for the documented methane/ethane state. The separate NRTL bubble/dew slice uses its own saved-source data and tolerance; neither gate is a general phase-envelope tracer.
- Phase-change heat capacity, general multiphase transport, and broad liquid-mixture transport are unsupported.

## Unit operations and flowsheets

- The public HTTP API exposes the versioned U0 routes only. Later unit-operation kernels are Python APIs until separate request schemas are approved.
- Flowsheet execution is deterministic for acyclic graphs. Material and energy recycles use the explicit bounded U8 solvers and are not automatically inferred from arbitrary graph cycles.
- Heat exchangers do not provide general phase-change rating, arbitrary shell layouts, mechanical design, or pressure/thermal co-iteration.
- Pipe thermal and hydraulic models support the documented saved segment conventions. Fully coupled transient multiphase networks, fitting catalogs, compressible networks, and relief discharge piping are unsupported.
- Reaction models cover the frozen stoichiometry, phase, kinetic basis, and thermodynamic closure named by each gate. General reaction-network selection, liquid Gibbs minimization, and fully coupled adiabatic kinetics are unsupported.
- Column functions include shortcut and frozen-profile material/energy gates plus simultaneous acetone/methanol NRTL total-condenser and reboiled-absorber solves. Both use a fixed pressure profile, one feed, a bottoms-flow specification, and calculated duties; the total condenser additionally uses a reflux-ratio specification. Other condenser/reboiler modes, pressure-profile generation, side draws, multiple feeds, alternate specifications, phase regimes, and arbitrary initial estimates are not a general rigorous-column package.

## Dynamics and controls

- Fixed-step explicit Euler is the deterministic regression integrator. Adaptive integration accepts only systems explicitly reduced to ODE form.
- `solve_ivp` is not presented as a general DAE solver. Models with unresolved algebraic states require an IDA-capable dependency that is not included in this release.
- Dynamic DWSIM parity is limited to the official water-tank/PID trajectory. General pressure-flow networks, multiphase vessel flashes, adaptive event location, actuator dynamics, and plant-wide initialization are unsupported.

## Runtime and deployment

- The service has no application authentication and is private-network only. Do not expose it directly to the public internet.
- Calculation requests use a two-process spawned worker pool. A timed-out calculation destroys that pool, which also aborts any other in-flight calculations assigned to it; the next request recreates the workers.
- The one-megabyte request limit and five-second calculation deadline are hard service boundaries.
- DWSIM is a Windows-only reference oracle and is not packaged with the Python service.
- Windows CPython 3.12.2 is locally verified. Linux amd64, Linux arm64, and macOS arm64 jobs are defined in the release workflow but must pass on the actual hosted runners before the 0.1.0 tag is created.

## Numerical expectations

- Tolerances are model- and reference-specific; they are recorded beside each golden in [compatibility.md](compatibility.md). A looser tolerance in one case must not be generalized to another.
- DWSIM's configured convergence criteria can leave a nonzero reference residual. The Python solver may close its own equations more tightly while using a documented parity tolerance for the saved reference.
- Results near phase boundaries, critical points, zero component flows, and correlation bounds require the domain checks exercised by the committed tests.
