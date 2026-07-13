# Lessons

- Measure generated Markdown fence counts before asserting an expected total.
- Choose and test the canonical SI base before writing molar or mass-flow factors.
- Treat entries from generic .NET dictionary enumerators as KeyValuePair-compatible objects, not DictionaryEntry.
- Verify the deployed script revision before diagnosing a repeated runtime error.
- Treat the script path in a stack trace as authoritative when verifying which artifact actually ran.
- When source text and a binding error conflict, capture the runtime invocation path and file hash before changing code.
- Do not prescribe an absolute clone path until the user has created that clone.
- Exclude capture timestamps from deterministic golden-case comparisons.
- Read Windows PowerShell JSON outputs with utf-8-sig to accept its UTF-8 BOM.
- Include -ExecutionPolicy Bypass whenever giving a direct PowerShell script invocation on the Windows capture host.
- Avoid PowerShell ETS property access on DWSIM objects because case-only CLR member duplicates break the adapter.
- Do not call GetType through a DWSIM PowerShell wrapper; unwrap BaseObject and invoke CLR members via System.Type.
- Read PSObject.BaseObject through reflected PSObject metadata because direct access is re-adapted to DWSIM under strict mode.
