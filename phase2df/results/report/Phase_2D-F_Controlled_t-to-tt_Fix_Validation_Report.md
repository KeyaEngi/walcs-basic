# Phase 2D-F Controlled t→tt Fix Validation Report

## 1. Validation scope

This was a controlled A/B validation of the two K2–K4 `t`→`tt` call-site corrections in `timmotion.f90`. The formal project was not used as a write target for the functional fix. SailModule, the adapter, `main.f90`, the wet-surface implementations, slamming, air-force, SectionLoad, and the project file were not functionally changed.

The four cases were:

| Case | Source | dt (s) | Simulated time (s) | Analysis window (s) |
|---|---|---:|---:|---:|
| B005 | baseline | 0.0050 | 0–7 | 3–7 |
| B0025 | baseline | 0.0025 | 0–7 | 3–7 |
| F005 | fixed | 0.0050 | 0–7 | 3–7 |
| F0025 | fixed | 0.0025 | 0–7 | 3–7 |

## 2. Baseline SHA-256

`baseline/timmotion.f90`:

`8A9C776BB08D9E2D5B9E00D307E9CC2FEDF918F939DA280E7C2F55F1382F4517`

This exactly matches the recovered formal-source reference hash.

## 3. Fixed SHA-256

`fixed/timmotion.f90`:

`2E4615127351E397A62335F0A35F27A1D8AA6B331CE9B723033822732BBFA42E`

The fixed file is two bytes longer than baseline, corresponding exactly to insertion of one additional `t` at each of the two permitted call sites.

## 4. Exact two-line diff

```diff
-            call Instant_wetsurface( y,it,t,smtf,Inst_ForceIS )
+            call Instant_wetsurface( y,it,tt,smtf,Inst_ForceIS )

-            call Instant_bdf_wetsurface( y,it,t,smtf,Instbdf_ForceIS )
+            call Instant_bdf_wetsurface( y,it,tt,smtf,Instbdf_ForceIS )
```

The K1 calls remain at full-step time `t`. Full-tree comparison of baseline and fixed found no other functional source difference.

## 5. Actual active nonlinear branch

The 605 input selects `Non_Linear = NL` and `nonlinearCtrl = 1`. The directly active wet-surface branch is therefore `Instant_bdf_wetsurface`.

K1 calls it with `t`; K2–K4 in the fixed copy call it with `tt = t + a(k)`.

## 6. Build status

The isolated fixed executable built and linked successfully with Intel Fortran Classic 2021.13, using the same Debug configuration and the same 59-object link set as the prior baseline build.

| Artifact | SHA-256 |
|---|---|
| Baseline executable | `F416786EC12214F8CF81F3EB874DF134C7DEC67BA1BE93F50F42757DBC33D8D4` |
| Fixed executable | `DD3C22BBB9006BFC2BDA5EADECDE40948E25AFEE71905DF5B0D64D22C4F00421` |

## 7. Run integrity

All four histories reach exactly 7.0 s, contain continuous Sail diagnostics, and have Sail enabled on every row.

| Case | History rows | Analysis rows | Final time (s) | Sail enabled |
|---|---:|---:|---:|---|
| B005 | 1400 | 801 | 7.0 | all rows |
| B0025 | 2800 | 1601 | 7.0 | all rows |
| F005 | 1400 | 801 | 7.0 | all rows |
| F0025 | 2800 | 1601 | 7.0 | all rows |

F005 and F0025 completed with exit code 0. No runtime error keywords were found.

## 8. Analysis method

The Phase 2C method was reused without alteration:

- analysis window: 3.0–7.0 s inclusive;
- encounter frequency: 2.52591142376 Hz;
- mean, population standard deviation, RMS, peak-to-peak;
- least-squares-equivalent single-frequency amplitude and phase after mean removal;
- phase comparisons use circular angular differences.

## 9. B005 result

| Variable | Mean | Std | RMS | P–P | Harmonic amp. | Phase (deg) |
|---|---:|---:|---:|---:|---:|---:|
| heave | -0.002224745 | 0.007549786 | 0.007870754 | 0.030284240 | 0.000400203 | 84.3208 |
| pitch | 0.038544782 | 0.000650424 | 0.038550269 | 0.002422787 | 0.000048246 | 47.7103 |
| V_REL_H_MAG | 33.576255594 | 0.008959128 | 33.576256789 | 0.043998928 | 0.000912353 | 82.8653 |
| Fx | -13280.595576 | 8.219161 | 13280.598119 | 40.365092 | 0.836880 | -97.1240 |
| My | -66402.977879 | 41.095807 | 66402.990596 | 201.825460 | 4.184398 | -97.1240 |

## 10. B0025 result

| Variable | Mean | Std | RMS | P–P | Harmonic amp. | Phase (deg) |
|---|---:|---:|---:|---:|---:|---:|
| heave | -0.002295102 | 0.006840083 | 0.007214861 | 0.029750358 | 0.000477830 | 81.6650 |
| pitch | 0.038355140 | 0.000517367 | 0.038358629 | 0.001936023 | 0.000032507 | 52.0964 |
| V_REL_H_MAG | 33.577028526 | 0.007910663 | 33.577029457 | 0.036555246 | 0.000551394 | 80.6210 |
| Fx | -13281.304422 | 7.257632 | 13281.306405 | 33.537085 | 0.505736 | -99.3667 |
| My | -66406.522110 | 36.288158 | 66406.532025 | 167.685423 | 2.528680 | -99.3667 |

## 11. F005 result

Every reported F005 value is identical to the corresponding B005 value at CSV precision. The entire files are byte-for-byte identical.

F005 history SHA-256:

`78218DBA85C0CFC2CE5E839F3469A211BB05594B332FAE67857EE823BEDB5AFE`

This is also the B005 history SHA-256.

## 12. F0025 result

Every reported F0025 value is identical to the corresponding B0025 value at CSV precision. The entire files are byte-for-byte identical.

F0025 history SHA-256:

`EA34D1B670A1851A2545272A654FCDFE7235A1516CCB3C599E4F731F5D7542D9`

This is also the B0025 history SHA-256.

## 13. Baseline convergence

Errors below are `|dt=.005 − dt=.0025| / |dt=.0025| × 100%`, except phase, which is an absolute circular difference in degrees.

| Variable | Mean | Std | RMS | P–P | Harmonic amp. | Phase |
|---|---:|---:|---:|---:|---:|---:|
| heave | 3.0656% | 10.3757% | 9.0909% | 1.7945% | 16.2457% | 2.6558° |
| pitch | 0.4944% | 25.7181% | 0.4996% | 25.1425% | 48.4164% | 4.3861° |
| V_REL_H_MAG | 0.002302% | 13.2538% | 0.002301% | 20.3628% | 65.4628% | 2.2442° |
| Fx | 0.005337% | 13.2485% | 0.005333% | 20.3596% | 65.4776% | 2.2427° |
| My | 0.005337% | 13.2485% | 0.005333% | 20.3596% | 65.4776% | 2.2427° |

## 14. Fixed convergence

The fixed convergence table is numerically identical to the baseline convergence table for every variable and every metric.

## 15. Convergence improvement

For all 30 variable/metric combinations:

`improvement = 0.0%`

Thus the K2–K4 `t`→`tt` correction produced neither amplitude improvement nor phase improvement in this controlled 605 experiment.

## 16. Pitch comparison

The two priority pitch measures are unchanged:

| Measure | Baseline dt difference | Fixed dt difference | Improvement |
|---|---:|---:|---:|
| pitch std | 25.7181% | 25.7181% | 0.0% |
| pitch harmonic amplitude | 48.4164% | 48.4164% | 0.0% |

The correction does not explain the remaining pitch-amplitude non-convergence.

## 17. V_REL comparison

| Measure | Baseline dt difference | Fixed dt difference | Improvement |
|---|---:|---:|---:|
| V_REL_H_MAG std | 13.2538% | 13.2538% | 0.0% |
| V_REL_H_MAG harmonic amplitude | 65.4628% | 65.4628% | 0.0% |

## 18. Fx comparison

| Measure | Baseline dt difference | Fixed dt difference | Improvement |
|---|---:|---:|---:|
| Fx std | 13.2485% | 13.2485% | 0.0% |
| Fx harmonic amplitude | 65.4776% | 65.4776% | 0.0% |

## 19. My comparison

| Measure | Baseline dt difference | Fixed dt difference | Improvement |
|---|---:|---:|---:|
| My std | 13.2485% | 13.2485% | 0.0% |
| My harmonic amplitude | 65.4776% | 65.4776% | 0.0% |

## 20. Phase comparison

No harmonic phase difference changed. The baseline/fixed coarse-to-fine phase differences remain:

- heave: 2.6558°;
- pitch: 4.3861°;
- V_REL_H_MAG: 2.2442°;
- Fx: 2.2427°;
- My: 2.2427°.

## 21. Why zero impact is still meaningful

`Instant_bdf_wetsurface` explicitly uses its time argument in regular- and irregular-wave phase expressions, so the original K2–K4 full-step-time call is a real stage time/state consistency defect. The fixed executable also has a different binary hash, proving the corrected source was rebuilt.

Nevertheless, the controlled output histories are byte-identical at both time steps. Therefore, for this active 605 path and these outputs, correcting that defect has no observable numerical effect. This experiment does not establish the internal reason for cancellation or insensitivity; it establishes the externally relevant A/B result.

## 22. Final causal assessment

The prescribed low-improvement conclusion applies:

> The t/tt mismatch is a confirmed algorithmic defect, but is not the dominant cause of the Phase 2C amplitude behavior.

It must not be described as the only cause or as a major contributor.

## 23. Merge recommendation

Recommend **MERGE** for the exact two-line correction because it restores RK stage time/state consistency and is narrowly scoped. This recommendation is based on algorithmic correctness, not on a convergence benefit: the controlled test demonstrated zero improvement in the current outputs.

Before merging, apply only the two reviewed lines to the formal tree and repeat its normal build/smoke checks. Do not copy unrelated isolated artifacts or backup files.

## 24. Phase 2B status

Phase 2B cannot be upgraded to PASS. The priority amplitude discrepancies remain unchanged, especially pitch harmonic amplitude (48.4164%) and the V_REL/Fx/My harmonic amplitudes (about 65.46–65.48%).

## 25. Formal parameter-study eligibility

The formal project is not yet eligible for a parameter study. This controlled fix removes one algorithmic defect but does not resolve or materially reduce Phase 2C amplitude non-convergence.

## 26. Execution note

An encoding-sensitive patch attempt briefly disturbed a working copy during preparation. It was immediately restored byte-for-byte from the intact source before the controlled build. End-state verification confirms the formal and baseline `timmotion.f90` hash is exactly the required reference hash, and the fixed tree contains only the permitted two functional line changes.

## 27. Final verdict

Phase 2D-F: **PASS**

t→tt formal fix recommendation: **MERGE**

Convergence impact: **NOT DOMINANT**

Phase 2B upgraded to PASS: **NO**

Formal parameter-study eligibility: **NO**
