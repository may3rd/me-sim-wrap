# Performance gate

Run the representative Phase 19 latency gate from the repository root:

```powershell
python benchmarks/benchmark_release.py
```

The script measures three direct kernel paths and the matching heavy HTTP worker boundaries. Its p95 limits are deliberately below half of the API's five-second calculation deadline. This is a release regression gate, not a substitute for load testing the deployed internal service.

Add native code only if a deployed workload misses an approved throughput target after Python algorithm and data-movement profiling. Keep the Python calculation as the numerical reference and run every native candidate against the same golden vectors.

The Phase 19 Windows reference run used CPython 3.12.2 with 20 core samples and five API samples. Observed p95 latencies were 3.22 ms for TP flash, 192.81 ms for a PH valve, 1.10 ms for the 121-point tank trajectory, 14.82 ms for the TP endpoint, 154.66 ms for the valve endpoint, and 332.93 ms for the U0 endpoint. These numbers are evidence for the current no-native-code decision, not portable performance guarantees.
