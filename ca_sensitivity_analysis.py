"""
CA sensitivity analysis for calibration.

Reads chronoamperometry (CA) traces exported as CSV (one file per electrode),
lets the user click the centre of each of the six current equilibriums
(Baseline + Drop1..Drop5), averages the current over a fixed 30 s window
centred on each click, and computes the sensor sensitivity from the linear
regression of the six average currents against the six known concentrations.

Output per run (written to the input folder):
  - sensitivity_analysis_<N>_files.csv   one row per electrode
  - sensor_response_<N>_files.pdf        the averaged plateaus, overlaid

Assumptions:
  - CSV columns 1 and 2 are time (s) and current (uA); extra columns ignored.
  - Data are acquired at 2 Hz  ->  a 30 s window is 60 samples.
  - Sensitivity = signed slope / electrode area  (uA mM^-1 cm^-2).
"""

import os
import glob
import csv

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# ---------------------------------------------------------------------------
# Configuration  (edit here if your setup differs)
# ---------------------------------------------------------------------------
CONCENTRATIONS_MM = [0.000, 0.019, 0.037, 0.054, 0.069, 0.083]  # x-values, 6 points
POINT_LABELS = ["Baseline", "Drop1", "Drop2", "Drop3", "Drop4", "Drop5"]
N_POINTS = len(CONCENTRATIONS_MM)               # 6

ELECTRODE_AREA_CM2 = 0.126
SAMPLING_RATE_HZ = 2                             # acquisition rate (fixed)
WINDOW_SECONDS = 30                             # averaging window per plateau
WINDOW_SAMPLES = WINDOW_SECONDS * SAMPLING_RATE_HZ   # 60 samples


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------
def get_folder_path():
    """Prompt for a folder until a valid directory is given."""
    while True:
        path = input("CSV folder path: ").strip().strip('"').strip("'")
        if os.path.isdir(path):
            return path
        print("  Not a valid folder, try again.")


def list_csv_files(folder):
    """Return the sorted list of .csv files in the folder."""
    return sorted(glob.glob(os.path.join(folder, "*.csv")))


def read_trace(filepath):
    """Read a CA trace; return (time, current) arrays, or (None, None) on error."""
    try:
        df = pd.read_csv(filepath).iloc[:, :2]          # first two columns only
        df.columns = ["time", "current"]
        df["time"] = pd.to_numeric(df["time"], errors="coerce")
        df["current"] = pd.to_numeric(df["current"], errors="coerce")
        df = df.dropna()
        return df["time"].values, df["current"].values
    except Exception as exc:
        print(f"  Could not read {os.path.basename(filepath)}: {exc}")
        return None, None


def clean_name(filepath):
    """File name without extension / trailing '_converted'."""
    name = os.path.splitext(os.path.basename(filepath))[0]
    if name.endswith("_converted"):
        name = name[: -len("_converted")]
    return name


# ---------------------------------------------------------------------------
# Plateau selection and averaging
# ---------------------------------------------------------------------------
def select_plateau_times(time, current, title):
    """Show the trace and let the user click the N_POINTS plateau centres.

    Returns the clicked x (time) values in click order. Extra clicks beyond
    N_POINTS are ignored.
    """
    picked = []

    fig, ax = plt.subplots(figsize=(15, 8))
    ax.plot(time, current, "b-", linewidth=1.5)
    ax.set_title(f"Click {N_POINTS} plateau centres (Baseline -> Drop5) - {title}")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Current (uA)")
    ax.grid(True, alpha=0.3)

    def on_click(event):
        if event.inaxes != ax or len(picked) >= N_POINTS:
            return
        picked.append(event.xdata)
        label = POINT_LABELS[len(picked) - 1]
        ax.axvline(event.xdata, color="red", ls="--", alpha=0.7)
        ax.text(event.xdata, ax.get_ylim()[1] * 0.95, label,
                ha="center", va="top", color="red", fontweight="bold", fontsize=8)
        fig.canvas.draw()
        print(f"  {label} at {event.xdata:.1f} s")
        if len(picked) == N_POINTS:
            print("  All points selected - close the window to continue.")

    cid = fig.canvas.mpl_connect("button_press_event", on_click)
    plt.show()
    fig.canvas.mpl_disconnect(cid)
    plt.close(fig)
    return picked


def window_around(time, current, center_time):
    """Mean current over a fixed 30 s window centred on center_time.

    Returns (window_current, mean_current). The window is WINDOW_SAMPLES
    points (60 at 2 Hz), clipped at the trace boundaries.
    """
    center_idx = int(np.argmin(np.abs(time - center_time)))
    half = WINDOW_SAMPLES // 2
    start = max(0, center_idx - half)
    end = min(len(current), center_idx + half)
    win = current[start:end]
    return win, float(np.mean(win))


# ---------------------------------------------------------------------------
# Per-file analysis
# ---------------------------------------------------------------------------
def analyse_file(filepath):
    """Process one CSV; return a results dict, or None if the file was skipped."""
    time, current = read_trace(filepath)
    if time is None:
        return None

    name = clean_name(filepath)
    print(f"\nFile: {name}")

    picked = select_plateau_times(time, current, name)
    if len(picked) != N_POINTS:
        print(f"  Skipped ({len(picked)} of {N_POINTS} points selected).")
        return None

    windows, means = [], []
    for label, t in zip(POINT_LABELS, picked):
        win, mean_i = window_around(time, current, t)
        windows.append(win)
        means.append(mean_i)
        print(f"  {label}: {mean_i:.3f} uA")

    # Linear fit: current (uA) vs concentration (mM)
    slope, intercept, r, p, _ = stats.linregress(CONCENTRATIONS_MM, means)
    sensitivity = slope / ELECTRODE_AREA_CM2          # signed, as requested

    return {
        "name": name,
        "windows": windows,
        "means": means,
        "sensitivity": sensitivity,
        "slope": slope,
        "intercept": intercept,
        "r2": r ** 2,
        "p": p,
    }


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
def save_plot(results, folder, n):
    """Overlay each electrode's six averaged plateaus and annotate concentrations."""
    plt.rcParams.update({"font.size": 12, "axes.linewidth": 1.2})
    fig, ax = plt.subplots(figsize=(12, 8), dpi=300)
    colors = plt.cm.tab10.colors

    for idx, res in enumerate(results):
        color = colors[idx % len(colors)]
        for j, win in enumerate(res["windows"]):
            # map each 30 s window onto a consecutive time slot: 0-30, 30-60, ...
            seg_t = np.linspace(j * WINDOW_SECONDS, (j + 1) * WINDOW_SECONDS, len(win))
            ax.plot(seg_t, win, color=color, linewidth=2.0, alpha=0.8,
                    label=res["name"] if j == 0 else None)

    y_top = ax.get_ylim()[1]
    for j, conc in enumerate(CONCENTRATIONS_MM):
        ax.text(j * WINDOW_SECONDS + WINDOW_SECONDS / 2, y_top,
                f"{conc:.3f} mM", ha="center", va="top", fontsize=10)

    ax.set_xlim(0, N_POINTS * WINDOW_SECONDS)
    ax.set_xticks(np.arange(0, N_POINTS * WINDOW_SECONDS + 1, WINDOW_SECONDS))
    ax.set_xlabel("Time (s)", fontweight="bold")
    ax.set_ylabel("Current (\u00b5A)", fontweight="bold")
    ax.set_title("Averaged plateaus (Baseline to Drop5)")
    ax.grid(True, alpha=0.2)
    ax.legend(loc="best", fontsize=11)
    fig.tight_layout()

    pdf_path = os.path.join(folder, f"sensor_response_{n}_files.pdf")
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\nPlot saved:   {pdf_path}")


def save_csv(results, folder, n):
    """Write one row per electrode to a summary CSV."""
    fields = (["File", "Sensitivity_uA_per_mM_cm2", "Slope_uA_per_mM",
               "Intercept_uA", "R_squared", "P_value"] +
              [f"{lbl}_uA" for lbl in POINT_LABELS])

    csv_path = os.path.join(folder, f"sensitivity_analysis_{n}_files.csv")
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for res in results:
            row = {
                "File": res["name"],
                "Sensitivity_uA_per_mM_cm2": res["sensitivity"],
                "Slope_uA_per_mM": res["slope"],
                "Intercept_uA": res["intercept"],
                "R_squared": res["r2"],
                "P_value": res["p"],
            }
            row.update({f"{lbl}_uA": v for lbl, v in zip(POINT_LABELS, res["means"])})
            writer.writerow(row)
    print(f"Results saved: {csv_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 60)
    print("CA sensitivity analysis")
    print("=" * 60)

    folder = get_folder_path()
    csv_files = list_csv_files(folder)
    if not csv_files:
        print("No CSV files found.")
        return
    print(f"Found {len(csv_files)} CSV file(s).")

    max_files = min(8, len(csv_files))
    while True:
        try:
            n = int(input(f"How many files to analyse (1-{max_files})? "))
            if 1 <= n <= max_files:
                break
        except ValueError:
            pass
        print(f"  Enter a number between 1 and {max_files}.")

    # Analyse each file; keep only the ones that were processed successfully.
    results = [r for r in (analyse_file(fp) for fp in csv_files[:n]) if r]
    if not results:
        print("\nNo files processed.")
        return

    print("\nSensitivity summary")
    print("-" * 60)
    for res in results:
        print(f"  {res['name']}: {res['sensitivity']:.3f} uA/mM/cm2, "
              f"R2 = {res['r2']:.4f}")

    save_plot(results, folder, len(results))
    save_csv(results, folder, len(results))
    print("\nDone.")


if __name__ == "__main__":
    main()
