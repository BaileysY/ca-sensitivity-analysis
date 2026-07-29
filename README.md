# CA Sensitivity Analysis

Custom script used to compute cholesterol-sensor sensitivity from
chronoamperometry (CA) calibration traces, as described in the
Supplementary Methods (Section S5).

## What it does

For each electrode (one CSV per electrode) the script:

1. Reads the time/current trace.
2. Lets you click the centre of the six current plateaus
   (Baseline, Drop1–Drop5).
3. Averages the current over a fixed **30 s window centred on each click**.
4. Linearly regresses (`scipy.stats.linregress`) the six average currents
   against the six known concentrations
   (0, 0.019, 0.037, 0.054, 0.069, 0.083 mM).
5. Reports **sensitivity = slope / electrode area (0.126 cm²)**.

## Requirements

- Python 3.8+
- numpy, pandas, scipy, matplotlib

```
pip install numpy pandas scipy matplotlib
```

## Input format

- One CSV per electrode, all in a single folder.
- Columns 1 and 2 = time (s) and current (µA); any extra columns are ignored.
- Data are assumed to be sampled at **2 Hz**, so 30 s = 60 points.

## Usage

```
python ca_sensitivity_analysis.py
```

Then follow the prompts:

1. Enter the folder path.
2. Enter how many CSVs to analyse (1–8).
3. For each file, click the six plateau centres in order
   (Baseline → Drop5), then close the plot window to continue.

## Output (written to the input folder)

- `sensitivity_analysis_<N>_files.csv` — one row per electrode:
  sensitivity, slope, intercept, R², p-value, and the six average currents.
- `sensor_response_<N>_files.pdf` — the averaged plateaus overlaid,
  annotated with concentration.

## Configuration

The parameters are grouped at the top of the script and can be edited if
your setup differs:

| Parameter            | Default                                   | Meaning                          |
|----------------------|-------------------------------------------|----------------------------------|
| `CONCENTRATIONS_MM`  | `[0, 0.019, 0.037, 0.054, 0.069, 0.083]`  | x-values for the regression (mM) |
| `ELECTRODE_AREA_CM2` | `0.126`                                   | geometric working-electrode area |
| `SAMPLING_RATE_HZ`   | `2`                                       | acquisition rate                 |
| `WINDOW_SECONDS`     | `30`                                      | averaging window per plateau     |

## Notes

- Sensitivity is the **signed** slope divided by the electrode area
  (µA·mM⁻¹·cm⁻²); take the absolute value if you only need the magnitude.
- The 30 s window is simply centred on each clicked point (no sub-window
  search).
