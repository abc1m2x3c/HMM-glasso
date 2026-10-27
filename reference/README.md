# Original author source

`src.R` is an unchanged copy of the source file in the author's **20250817 mod 8 corrected** simulation directory. It contains `Graph.gen.multiK` (including chain/mod 5 and random/mod 8), `EM.glasso`, `hierarchical_gmm3`, and earlier experimental functions.

SHA-256:

```
0A3FD1CD1DB9FD974DB3E3BA0CD1FD3F328B09FD5587732C0328B55E02982679
```

The corresponding **20250803 mod 5 corrected** source has the same flat-mixture and hierarchical fitting statements, with additional iteration logging. The portable demo is based on these common updates.

This reference file is archived for provenance and is **not sourced by the demo**. Some functions require unqualified packages or global tuning variables. It also contains historical evaluation and BIC code; use the explicitly labeled routines in `R/` for the showcase. See [implementation notes](../docs/implementation.md) for the refactor and interpretation details.
