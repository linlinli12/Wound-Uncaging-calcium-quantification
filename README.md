# Image analysis of wound- and ATP-uncaging-induced Ca²⁺ waves

Code used to measure the Ca²⁺ waves in the fin epithelium for the paper: the wave triggered by a laser wound and
the wave triggered by ATP uncaging. Both are measured with the same method:

- **Sheet pixels.** Besides the wave, single cells flash on their own before, during and after it. Only the pixels
  reached by the sheet of responding cells that grows out from the source (the wound edge or the lit uncaging
  circle) are measured; cells that flash on their own are left out (`ca_wave/sheet_wave.py`, shared by both
  analyses).
- **Sectors.** Each movie is measured in 11.25° wedges × 25-µm rings around the wound or the uncaging spot, over the
  full circle. Each sector's trace (the mean of its sheet pixels) is normalized as `(I - min)/min`, and its arrival
  time is the half-max crossing. A sector counts as a response only if its trace rises clearly above its own noise.
- **Wave test.** An outward wave is reported only if the responses reach far enough, travel at least 50 µm beyond
  the immediate response (the rings that respond within 21 s), and arrive later with distance (ring-label
  permutation test). The wave speed is a Theil-Sen fit of the median arrival time per ring against distance.

## Pipeline

| Step | File | Input | Output |
|---|---|---|---|
| 1 | `notebooks/1_wound_hole_tracking.ipynb` | membrane TIFFs of one wound dataset | `wound_tracking/<membrane folder>/hole_tracking_results_new.csv` (wound centre and radius per frame) |
| 2 | `notebooks/2_wound_wave_analysis.ipynb` | calcium + membrane TIFFs and the wound tracking of one dataset | `results/wound_wave_analysis/<dataset>/` |
| 3 | `notebooks/3_uncaging_wave_analysis.ipynb` | uncaging AVIs (all stacks) | `results/uncaging_wave_analysis/` |
| 4 | `scripts/4_data_timing.py` | outputs of 2 (every wound dataset) and 3, uncaging AVIs | `results/data_timing_by_ring.csv` (onset, half-max and peak time per ring) |
| 5 | `scripts/5_export_group_targets.py` | outputs of 2-4 | `results/group_targets.csv` (group means per ring, used to fit the model) |

Notebooks 2 and 3 import the sheet-pixel code from `ca_wave/sheet_wave.py` (the method settings used in the paper
are its `SHEET_SETTINGS`); run them from the `notebooks/` folder so that the repository root is found.

Notebooks 1 and 2 run on **one wound dataset** at a time; they are set up for the example dataset (below) and each
explains its method step by step. To analyse another wound dataset, set its folders and settings in both notebooks
(notebook 2 lists the settings used for every dataset in the paper). Steps 4 and 5 combine all datasets, so they
need the full image data, with notebook 2 run once for every wound dataset.

## Acquisition and time calibration

- **1 px = 1 µm.**
- **One frame = 7 s**, for both wound and uncaging movies. Each analysed frame is the average of 10 raw frames of
  0.7 s; an uncaging movie has 40 frames (400 raw frames, 280 s).
- The TIFFs carry no time calibration, and the frame rate in the AVI headers (10 fps) is **not** the acquisition
  rate. All code reads frames by index and converts them with `FRAME_INTERVAL_S = 7.0`.
- Uncaging: the light is on from frame 0; in the "35 s" stacks it goes off after raw frame 50 (AVG frame 5).
- Wound: time 0 is the wound frame, the frame in which the membrane inside the hole first turns dark (set per
  dataset as `wound_frame` in notebook 2).

## Example data

One wound dataset, **ctrl_0820** (control), is included in `data/wound/` (65 frames per channel, about 130 MB).
Notebooks 1 and 2 are set up for it and run as they are, in this order:

```
pip install -r requirements.txt
cd notebooks
jupyter nbconvert --to notebook --execute --inplace 1_wound_hole_tracking.ipynb
jupyter nbconvert --to notebook --execute --inplace 2_wound_wave_analysis.ipynb
```

- **Notebook 1** tracks the wound with the settings used for the paper (`METHOD = "halfway"`) and writes
  `wound_tracking/ctrl_0820_membrane/hole_tracking_results_new.csv`: the wound radius goes from 40 px at frame 21
  to 10 px at frame 64.
- **Notebook 2** reads that tracking, measures the Ca²⁺ wave and writes its tables and figures to
  `results/wound_wave_analysis/ctrl_0820/`. Expected result: an outward wave reaching 200 µm, immediate response
  25 µm, Kendall τ = 0.56 (p = 0.0002), speed 2.43 µm/s (95 % CI 1.39-4.09).

## Using your own data

Each notebook starts with a **"Settings to edit"** cell and a **checklist for your own data**; the method settings
used in the paper are kept in a separate cell, so they stay unchanged unless you choose to change them.

- **Wound movies** (notebooks 1 and 2): for each movie, a folder with the Ca²⁺ channel and a folder with a
  channel in which the wound is a dark hole (e.g. a membrane marker), one 2-D TIFF per frame, pixel-aligned and
  with the same frames. Track the wound with notebook 1 (set `MEMBRANE_FOLDER`, `START_FRAME` and the seed box
  `INIT_ROI`; the notebook draws the box on the image), then set `DATASET` in notebook 2 (folders, `wound_frame`,
  `angles`) and run it. Each notebook has a troubleshooting table.
- **Uncaging movies** (notebook 3): one grayscale AVI per experiment; add an entry to `STACKS`.
- **Pixel size and frame interval**: set `PIXEL_SIZE_UM` and `FRAME_INTERVAL_S` for your microscope. The rings are
  defined in pixels (25 px apart, out to 200 px); if your pixel size is not 1 µm, set the ring sizes as the
  checklists describe to keep 25-µm rings.

Keep the method settings of the paper to compare your results with the paper's.

## Data layout

Apart from the example dataset, the image data are not part of this repository (see the paper's Data
Availability statement). Put them in `data/` as below, or point the environment variables to where they are.

```
data/
  uncaging/                                   # UNCAGING_DATA_DIR
    101626 ATP uncaging (AVG_Stack)/
      uncage time (full)/AVG_Stack_ATP1_2.avi
      uncage time (35s or frame 50-400)/AVG_Stack_ATP3.avi, AVG_Stack_ATP4.avi, AVG_Stack_ATP5.avi
    100126 uncaging suramin (AVG_Stack)/
      AVG_Stack_Suramin1.avi, AVG_Stack_ATP uncaging suramin 2.avi, AVG_Stack_ATP uncaging suramin 3.avi
  wound/                                      # WOUND_DATA_DIR
    ctrl_0808_calcium/AVG_Stack0000.tif ...   ctrl_0808_membrane/MAX_Stack0000.tif ...
    ctrl_0820_calcium/   ctrl_0820_membrane/
    pom1-3_calcium/      pom1-3_membrane/
    Gd3_3_calcium/       Gd3_3_membrane/
    Gd3-4-calcium/       Gd3-4-membrane/
    CBX_5_calcium/       CBX_5_membrane/        (likewise CBX_6, CBX_7)
    hep_6_calcium/       hep_6_membrane/        (likewise hep_7)
    heptonal_3_calcium/  heptonal_3_membrane/
    heptonal_5_calcium/  heptanol_5_membrane/
```

Each wound folder holds one TIFF per frame (calcium: `AVG_Stack####.tif`; membrane: `AVG_Stack####.tif` or
`MAX_Stack####.tif`); the two channels of a dataset are pixel-aligned and their frames match one to one. If an
uncaging AVI is not at its listed path, notebook 3 searches for it by file name anywhere under `data/uncaging/`.

| Environment variable | Default | Used by |
|---|---|---|
| `UNCAGING_DATA_DIR` | `data/uncaging` | notebook 3, script 4 |
| `WOUND_DATA_DIR` | `data/wound` | notebooks 1 and 2 |
| `WOUND_TRACKING_DIR` | `wound_tracking` | notebooks 1 and 2 |
| `RESULTS_DIR` | `results` | notebooks 2 and 3, scripts 4 and 5 |

## How to run the full analysis

With all image data in `data/` (see "Data layout"):

1. For every wound dataset, run notebook 2 with its settings (table in notebook 2), using the supplied tracking
   (`TRACKING_CSV = "hole_tracking_results.csv"`) or, for the five `"halfway"` datasets, the tracking from
   notebook 1 run with that dataset's settings (table in notebook 1).
2. Run notebook 3 (all uncaging stacks at once).
3. Run the two scripts:

```
python3 scripts/4_data_timing.py
python3 scripts/5_export_group_targets.py group_targets.csv --source --onset
```

Open the notebooks in Jupyter and run all cells from the `notebooks/` folder, or run them with
`jupyter nbconvert --execute` as above. Each notebook takes about a minute.

## Wound tracking

Notebook 1 finds the wound (a dark hole in the membrane channel) on a start frame and follows it frame by frame as
it closes. It has two methods (`METHOD`):

- **`"halfway"`**: the hole edge is put where the membrane is halfway between the hole floor and the surrounding
  tissue, both measured again in every frame, so the edge follows the hole as it fills in and closes. Used for
  **ctrl_0808, ctrl_0820, Gd3-4, CBX_6 and hep_6**; their settings are listed in notebook 1, and notebook 1
  reproduces their tracking exactly.
- **`"otsu"`**: the hole is the largest region below an Otsu threshold inside a search window that follows the
  hole. Used for the **other wound datasets**. Once a filling hole is only a little darker than the tissue, this
  region spreads beyond the hole and the radius is held; that is why the five datasets above were re-tracked.

Image quality differs between datasets, so the settings of the `"otsu"` method (`START_FRAME`, `INIT_ROI`,
`SWITCH_FRAME`, window sizes) were tuned for each dataset. **These per-dataset settings were not recorded.** The
resulting tracking of every dataset is supplied in `wound_tracking/<membrane folder>/hole_tracking_results.csv`.
Notebook 1 writes `hole_tracking_results_new.csv` and never overwrites the supplied file.

## Datasets

| Experiment | Group | Datasets |
|---|---|---|
| Wound | control | ctrl_0808, ctrl_0820 |
| | POM1 | pom1-3 |
| | Gd³⁺ | Gd3_3, Gd3-4 |
| | carbenoxolone (CBX) | CBX_5, CBX_6, CBX_7 |
| | heptanol | hep_6, hep_7, heptanol_3, heptanol_5 |
| ATP uncaging | ATP, light on for the whole movie (280 s) | ATP1_2 |
| | ATP, light on for 35 s | ATP3, ATP4, ATP5 |
| | ATP + suramin, light on for the whole movie | Suramin1, suramin 2, suramin 3 |

The wound dataset suramin_1 is excluded from the analysis (data-quality issue).

## Requirements

Python ≥ 3.10 and the packages in `requirements.txt`. Tested with numpy 2.2, scipy 1.15, pandas 2.3,
matplotlib 3.10, scikit-image 0.25, opencv-python 5.0 and tifffile 2025.5.

## Citation and license

Please cite the paper when using this code: *[citation to be added]*.
License: *[to be added]*.
