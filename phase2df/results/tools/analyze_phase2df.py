from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"D:\desktop\mxf毕业论文\walcas 源代码\phase2df\results")
OUT = ROOT / "report"
CASES = {
    "B005": Path(r"C:\Users\34675\.codex\visualizations\2026\08\12\019ff3f4-ba37-7b82-be6c-2d142765800b\phase2c_results\dt_005\tempout\605_LC0001_sail_history.csv"),
    "B0025": Path(r"C:\Users\34675\.codex\visualizations\2026\08\12\019ff3f4-ba37-7b82-be6c-2d142765800b\phase2c_results\dt_0025\tempout\605_LC0001_sail_history.csv"),
    "F005": Path(r"C:\Users\34675\.codex\visualizations\2026\08\13\019ff913-e210-7f62-97d8-5da069f46b83\phase2df_runs\F005\tempout\605_LC0001_sail_history.csv"),
    "F0025": Path(r"C:\Users\34675\.codex\visualizations\2026\08\13\019ff913-e210-7f62-97d8-5da069f46b83\phase2df_runs\F0025\tempout\605_LC0001_sail_history.csv"),
}
VARIABLES = {
    "heave": "heave",
    "pitch": "pitch",
    "V_REL_H_MAG": "V_REL_H_MAG",
    "Fx": "FX_BODY",
    "My": "MY_BODY",
}
METRICS = ("mean", "std", "rms", "ptp", "harmonic_amplitude", "harmonic_phase_deg")
U0 = 32.07798433303833
OMEGA = 2.36
G = 9.8
K = OMEGA**2 / G
OMEGA_E = OMEGA - K * U0
F_E = abs(OMEGA_E) / (2.0 * math.pi)
T0 = 3.0
T1 = 7.0


def circular_difference_deg(a: float, b: float) -> float:
    return abs((a - b + 180.0) % 360.0 - 180.0)


def stats(t: np.ndarray, x: np.ndarray) -> dict[str, float]:
    xc = x - np.mean(x)
    angle = 2.0 * math.pi * F_E * t
    aa = 2.0 / len(x) * np.sum(xc * np.cos(angle))
    bb = 2.0 / len(x) * np.sum(xc * np.sin(angle))
    return {
        "mean": float(np.mean(x)),
        "std": float(np.std(x, ddof=0)),
        "rms": float(np.sqrt(np.mean(x * x))),
        "ptp": float(np.ptp(x)),
        "harmonic_amplitude": float(math.hypot(aa, bb)),
        "harmonic_phase_deg": float(math.degrees(math.atan2(-bb, aa))),
    }


def relative_difference(coarse: float, fine: float) -> float:
    scale = max(abs(fine), 1.0e-30)
    return abs(coarse - fine) / scale * 100.0


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    integrity: list[dict[str, object]] = []
    all_stats: dict[str, dict[str, dict[str, float]]] = {}

    for case, path in CASES.items():
        df = pd.read_csv(path)
        t_all = df["time"].to_numpy(dtype=float)
        mask = (t_all >= T0 - 1.0e-12) & (t_all <= T1 + 1.0e-12)
        t = t_all[mask]
        all_stats[case] = {}
        integrity.append({
            "case": case,
            "rows": len(df),
            "time_first": float(t_all[0]),
            "time_last": float(t_all[-1]),
            "window_rows": int(mask.sum()),
            "sail_all_enabled": bool(df["sail_enable"].astype(str).str.upper().eq("T").all()),
        })
        for variable, column in VARIABLES.items():
            result = stats(t, df.loc[mask, column].to_numpy(dtype=float))
            all_stats[case][variable] = result
            for metric, value in result.items():
                rows.append({"case": case, "variable": variable, "metric": metric, "value": value})

    convergence: list[dict[str, object]] = []
    for variable in VARIABLES:
        for metric in METRICS:
            b0 = all_stats["B005"][variable][metric]
            b1 = all_stats["B0025"][variable][metric]
            f0 = all_stats["F005"][variable][metric]
            f1 = all_stats["F0025"][variable][metric]
            if metric == "harmonic_phase_deg":
                be = circular_difference_deg(b0, b1)
                fe = circular_difference_deg(f0, f1)
                unit = "deg"
            else:
                be = relative_difference(b0, b1)
                fe = relative_difference(f0, f1)
                unit = "percent"
            improvement = (be - fe) / be * 100.0 if be > 1.0e-30 else (0.0 if fe <= 1.0e-30 else -math.inf)
            convergence.append({
                "variable": variable,
                "metric": metric,
                "baseline_error": be,
                "fixed_error": fe,
                "error_unit": unit,
                "improvement_percent": improvement,
            })

    pd.DataFrame(rows).to_csv(OUT / "phase2df_metrics.csv", index=False)
    pd.DataFrame(integrity).to_csv(OUT / "phase2df_integrity.csv", index=False)
    pd.DataFrame(convergence).to_csv(OUT / "phase2df_convergence.csv", index=False)

    print(f"encounter_frequency_hz={F_E:.12g}")
    print("INTEGRITY")
    print(pd.DataFrame(integrity).to_string(index=False))
    print("METRICS")
    for variable in VARIABLES:
        print(f"[{variable}]")
        print("case,mean,std,rms,ptp,harmonic_amplitude,harmonic_phase_deg")
        for case in CASES:
            s = all_stats[case][variable]
            print(case + "," + ",".join(f"{s[m]:.12g}" for m in METRICS))
    print("CONVERGENCE")
    print(pd.DataFrame(convergence).to_csv(index=False, quoting=csv.QUOTE_MINIMAL))


if __name__ == "__main__":
    main()
