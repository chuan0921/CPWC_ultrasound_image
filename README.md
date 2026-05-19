# CPWC Ultrasound Image Simulation

This project generates synthetic CPWC-style ultrasound image sequences for
speckle tracking and WSS/WSR algorithm testing.

The current primary workflow is the fast image-domain simulator in
`main_flow_simulation.m`. It does not run Field II. Field II code is kept only
for small sanity checks.

## Main Workflow

Run the fast simulator from MATLAB:

```matlab
main_flow_simulation
```

The default probe profile is `literature_l12_3v`. To run another probe profile,
set `probe_name` before calling the simulator:

```matlab
probe_name = 'literature_l12_3v';
main_flow_simulation

probe_name = 'linear_l12_5mhz';
main_flow_simulation

probe_name = 'zipper_array';
main_flow_simulation
```

Each run creates a timestamped output folder:

```text
output/<probe_name>/run_yyyymmdd_HHMMSS_mmm/
```

For the default settings this is:

```text
output/literature_l12_3v/run_.../
```

Each seed is saved separately:

```text
seed_0001/
seed_0002/
...
seed_0005/
```

## Tracking Input

For linear-array profiles, the official tracking input is:

```matlab
lri_env_frames(:,:,angle_idx,frame_idx)
```

The data layout is fixed as:

```text
[z, x, angle, frame]
```

For the default literature probe settings:

```matlab
size(lri_env_frames)
% [167 104 5 20]
```

Use the same angle across consecutive frames for frame-to-frame tracking:

```matlab
load('image_data.mat', 'lri_env_frames', 'params')

angle_idx = 3;  % 0 deg
img1 = lri_env_frames(:,:,angle_idx,1);
img2 = lri_env_frames(:,:,angle_idx,2);
```

Timing and pixel size:

```matlab
dt_frame = params.n_angles / params.PRF;
dx = params.dx;
dz = params.dz;
```

Do not use compounded images for tracking. This project intentionally saves
uncompounded per-angle LRI envelope images.

For the zipper-array profile, the official row-specific tracking input is:

```matlab
lri_env_rows(:,:,row_idx,angle_idx,frame_idx)
```

The zipper data layout is:

```text
[z, x, row, angle, frame]
```

Use `row_idx = 1` and `row_idx = 2` to access the two elevation rows. The
combined `lri_env_frames` output is kept only for preview/backward-compatible
inspection and should not be used for zipper row-to-row tracking.

## Seed Folder Contents

Each `seed_XXXX` folder contains only tracking-related files:

```text
image_data.mat
ground_truth.mat
tracking_metadata.mat
lri_frames/
lri_angles_frame1.png
```

`image_data.mat` contains:

```text
lri_env_frames
lri_bmode_frames
params
timing
tracking_metadata
```

For zipper-array runs, `image_data.mat` also contains:

```text
lri_env_rows
lri_bmode_rows
```

`lri_frames/` contains per-frame/per-angle `.mat` and `.png` files. The `.mat`
files store one LRI envelope image and its B-mode visualization:

```text
frame_001_angle_+00deg.mat
frame_001_angle_+00deg.png
...
```

The `.png` files are for visual inspection only. Quantitative tracking should
use `lri_env_frames` from `image_data.mat`.

## Default Simulation Settings

The default settings are defined in `setup_parameters.m`.

Key defaults:

```text
Probe: literature_l12_3v
Center frequency: 8 MHz
Angles: [-10 -5 0 5 10] deg
PRF: 10 kHz
Seeds: 1:5
Frames per seed: 20
Velocity profile: steady Poiseuille flow
Center velocity: 0.50 m/s
Vessel radius: 2.5 mm
Vessel inner diameter: 5 mm
Vessel outer diameter: 5.4 mm
Vessel center depth: 20 mm
Wall thickness: 0.2 mm
Elevation support: +/-5 mm for linear-array profiles
Elevation beam sigma: probe-derived from wavelength, elevation focus, and element height
Noise mode: SNR-scaled complex image noise
```

The vessel radius defines the lumen. The wall is added outside the lumen.
For linear-array profiles, the elevation simulation support and elevation
projection weighting are fixed by probe-related settings rather than vessel
radius, so changing the lumen radius does not also change the elevation beam
model.

## Ground Truth

Ground truth is saved in `ground_truth.mat`.

Important fields:

```matlab
ground_truth.flow.vx_map
ground_truth.flow.vessel_mask
ground_truth.flow.dx_map
ground_truth.flow.dx_max
ground_truth.vessel.center_z
ground_truth.vessel.radius
ground_truth.vessel.wall_upper
ground_truth.vessel.wall_lower
```

For the default setup:

```text
vessel center = 20 mm
upper lumen boundary = 17.5 mm
lower lumen boundary = 22.5 mm
upper outer wall boundary = 17.3 mm
lower outer wall boundary = 22.7 mm
maximum frame displacement = 0.25 mm/frame
```

## Simulation Model

The current fast simulator is image-domain based:

1. Generate filtered complex speckle fields for tissue, flow, and wall.
2. Apply the analytical Poiseuille displacement field to the flow speckle.
3. Generate one uncompounded LRI envelope image per angle and frame.
4. Save only per-angle LRI envelope data for tracking.

This is not a full Field II RF/channel simulation. It is intended for fast
algorithm development and controlled WSS/WSR testing.

## Field II Sanity Check

`main_fieldii_parallel_simulation.m` and `simulate_fieldii_seed.m` are retained
for small Field II sanity checks with the literature linear array.

The Field II sanity-check script reads the local Field II MATLAB folder from
the `FIELDII_DIR` environment variable. Do not commit local absolute paths.

Example:

```bash
export FIELDII_DIR=/path/to/fieldii/m_files
```

The Field II output uses the same tracking interface:

```matlab
lri_env_frames(:,:,angle_idx,frame_idx)
```

## Probe Profiles

Probe settings are centralized in `setup_parameters.m`.

Available profiles:

```text
literature_l12_3v   Literature L12-3v, 8 MHz, 128-element linear array
linear_l12_5mhz     L12 geometry, 5 MHz, 128-element linear array
zipper_array        5 MHz zipper array, 128 x 2 = 256 elements
```

To select a profile, set `probe_name` before running `main_flow_simulation.m`.
The simulator passes that name to `setup_parameters(probe_name)`.

The fast simulator can use these settings for image-domain simulation. The
Field II sanity check currently supports only the linear array profile.

### Adding a New Array

Add new probe profiles in `setup_parameters.m`, inside the `probe_profile(name)`
function. Each profile should define a new `case` with the required probe fields:

```matlab
probe.name
probe.label
probe.geometry
probe.f0
probe.n_elements
probe.pitch
probe.kerf
probe.element_width
probe.element_height
probe.elev_focus
probe.n_sub_x
probe.n_sub_y
probe.fs
```

For a standard linear array, these scalar geometry fields are enough. For a
non-standard geometry such as the zipper array, also define:

```matlab
probe.element_data
probe.element_centers
```

`element_data` uses the Field-II-style element polygon format already used by
the zipper-array profile. After adding the profile, run it by setting:

```matlab
probe_name = '<new_profile_name>';
main_flow_simulation
```

## Version Control

Generated data are ignored by git:

```text
output/
output_fieldii_parallel/
```

Reference PDFs and local editor settings are also ignored.
