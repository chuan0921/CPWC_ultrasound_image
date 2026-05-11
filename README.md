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

The official tracking input is:

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
Vessel radius: 3 mm
Vessel center depth: 20 mm
Wall thickness: 0.2 mm
```

The vessel radius defines the lumen. The wall is added outside the lumen.

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
upper lumen boundary = 17 mm
lower lumen boundary = 23 mm
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

The Field II path is currently set in `main_fieldii_parallel_simulation.m`:

```matlab
fieldii_dir = '/Users/dorishsu/MBP16/NTU/PCLAB/Project/m_files';
```

The Field II output uses the same tracking interface:

```matlab
lri_env_frames(:,:,angle_idx,frame_idx)
```

## Probe Profiles

Probe settings are centralized in `setup_parameters.m`.

Available profiles:

```matlab
setup_parameters('literature_l12_3v')
setup_parameters('zipper_array')
```

The fast simulator can use these settings for image-domain simulation. The
Field II sanity check currently supports only the linear array profile.

## Version Control

Generated data are ignored by git:

```text
output/
output_fieldii_parallel/
```

Reference PDFs and local editor settings are also ignored.
