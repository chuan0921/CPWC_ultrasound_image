%% main_flow_simulation.m
% Fast image-domain CEUS flow phantom.
% Goal: reproduce the simulated B-mode appearance without Field II RF.

clear; close all; clc;

params = setup_parameters();

proj_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(proj_dir, 'output');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
probe_out_dir = fullfile(out_dir, params.probe.name);
if ~exist(probe_out_dir, 'dir'), mkdir(probe_out_dir); end
run_id = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
run_out_dir = fullfile(probe_out_dir, ['run_' run_id]);
if ~exist(run_out_dir, 'dir'), mkdir(run_out_dir); end
params.run_id = run_id;
params.run_timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
params.output_root = run_out_dir;

fprintf('=== Fast CEUS Image Simulation Batch ===\n');
fprintf('Probe: %s (%s)\n', params.probe.name, params.probe.label);
fprintf('Grid: %d x %d pixels, frames/seed=%d, seeds=%d\n', ...
    params.Nz, params.Nx, params.n_frames, numel(params.random_seeds));
fprintf('Run ID: %s\n', params.run_id);
fprintf('Output: %s\n', run_out_dir);

seed_results = cell(numel(params.random_seeds), 1);

for seed_idx = 1:numel(params.random_seeds)
    params.random_seed = params.random_seeds(seed_idx);
    rng(params.random_seed);

    seed_dir = fullfile(run_out_dir, sprintf('seed_%04d', params.random_seed));
    if ~exist(seed_dir, 'dir'), mkdir(seed_dir); end

    fprintf('--- Seed %d (%d/%d) ---\n', ...
        params.random_seed, seed_idx, numel(params.random_seeds));

    tic_seed = tic;
    [lri_env_frames, lri_bmode_frames] = simulate_ceus_images(params);
    timing.total = toc(tic_seed);
    timing.seed = params.random_seed;
    fprintf('Image generation: %.2f s\n', timing.total);

    ground_truth = scatterer_medium('ground_truth', params);
    tracking_metadata = build_tracking_metadata(params, 'fast_image_domain');

    save(fullfile(seed_dir, 'image_data.mat'), ...
        'lri_env_frames', 'lri_bmode_frames', ...
        'params', 'timing', 'tracking_metadata', '-v7.3');
    save(fullfile(seed_dir, 'ground_truth.mat'), 'ground_truth', '-v7.3');
    save(fullfile(seed_dir, 'tracking_metadata.mat'), 'tracking_metadata', '-v7.3');

    save_tracking_images(params, lri_env_frames, lri_bmode_frames, tracking_metadata, seed_dir);
    save_lri_montage(params, lri_bmode_frames, seed_dir);
    close all;

    seed_results{seed_idx}.seed = params.random_seed;
    seed_results{seed_idx}.out_dir = seed_dir;
    seed_results{seed_idx}.total_seconds = timing.total;
    fprintf('Data saved to: %s\n', seed_dir);
end

save(fullfile(run_out_dir, 'batch_summary.mat'), 'seed_results', 'params', '-v7.3');
fprintf('=== Done ===\n');

function [lri_env_frames, lri_bmode_frames] = simulate_ceus_images(params)
    [X, Z] = meshgrid(params.x_grid, params.z_grid);

    z0 = params.vessel_center_z;
    r = abs(Z - z0);
    lumen = r <= params.R;
    wall = r > params.R & r <= params.R + params.wall_thickness;
    tissue = ~lumen & ~wall;

    % Image-domain point spread function. This approximates compounded plane
    % wave resolution without simulating channel RF.
    axial_psf = gaussian_kernel(params.lambda / 3, params.dz);
    lateral_psf = gaussian_kernel(params.lambda * 0.85, params.dx);

    tissue_iq = filtered_speckle(params, axial_psf, lateral_psf);
    flow_iq0 = filtered_speckle(params, axial_psf, lateral_psf);
    wall_iq = filtered_speckle(params, axial_psf, lateral_psf);

    depth_gain = exp(-14 * (Z - params.z_min));
    lateral_gain = exp(-0.5 * (X / (0.42 * max(abs(params.x_grid)))).^2);
    sensitivity = depth_gain .* (0.55 + 0.45 * lateral_gain);

    tissue_iq = 0.28 * tissue_iq .* tissue;
    flow_iq0 = 1.25 * flow_iq0 .* lumen;
    wall_iq = 2.6 * wall_iq .* wall;

    vx_map = zeros(params.Nz, params.Nx);
    vx_map(lumen) = params.v0 * (1 - (r(lumen) / params.R).^2);

    lri_env_frames = zeros(params.Nz, params.Nx, params.n_angles, params.n_frames);
    lri_bmode_frames = zeros(params.Nz, params.Nx, params.n_angles, params.n_frames);

    for frame = 1:params.n_frames
        for angle_idx = 1:params.n_angles
            pulse_idx = (frame - 1) * params.n_angles + (angle_idx - 1);
            displacement = vx_map * pulse_idx / params.PRF;
            flow_iq = interp2(X, Z, flow_iq0, X - displacement, Z, 'linear', 0);

            angle = params.angles(angle_idx);
            max_angle = max(abs(params.angles));
            if max_angle > 0
                angle_gain = 1 - 0.18 * abs(angle) / max_angle;
            else
                angle_gain = 1;
            end
            angle_noise = params.image_noise_floor * complex(randn(size(X)), randn(size(X)));

            shadow = 1 - 0.18 * lumen;
            lri_iq = angle_gain * (tissue_iq + flow_iq + wall_iq) .* sensitivity .* shadow;
            lri_iq = lri_iq + angle_noise;

            lri_env_frames(:,:,angle_idx,frame) = abs(lri_iq);
            lri_bmode_frames(:,:,angle_idx,frame) = ...
                envelope_to_bmode(lri_env_frames(:,:,angle_idx,frame), params.dynamic_range);
        end
    end

    expected_size = [params.Nz, params.Nx, params.n_angles, params.n_frames];
    assert(isequal(size(lri_env_frames), expected_size), ...
        'lri_env_frames must be [z, x, angle, frame].');
end

function iq = filtered_speckle(params, axial_psf, lateral_psf)
    iq = complex(randn(params.Nz, params.Nx), randn(params.Nz, params.Nx));
    iq = conv2(conv2(iq, axial_psf(:), 'same'), lateral_psf(:).', 'same');
end

function kernel = gaussian_kernel(sigma_m, spacing_m)
    radius = max(1, ceil(4 * sigma_m / spacing_m));
    x = (-radius:radius) * spacing_m;
    kernel = exp(-0.5 * (x / sigma_m).^2);
    kernel = kernel / sum(kernel);
end

function bmode = envelope_to_bmode(envelope, dynamic_range)
    env_max = max(envelope(:));
    bmode = 20 * log10(max(envelope, eps) / max(env_max, eps));
    bmode = max(bmode, -dynamic_range);
end

function save_tracking_images(params, lri_env_frames, lri_bmode_frames, metadata, out_dir)
    lri_dir = fullfile(out_dir, 'lri_frames');
    if ~exist(lri_dir, 'dir'), mkdir(lri_dir); end

    for frame = 1:params.n_frames
        for angle_idx = 1:params.n_angles
            angle_deg = round(rad2deg(params.angles(angle_idx)));
            env_img = lri_env_frames(:,:,angle_idx,frame);
            bmode_img = lri_bmode_frames(:,:,angle_idx,frame);
            png_img = bmode_to_uint8(bmode_img, params.dynamic_range);

            stem = sprintf('frame_%03d_angle_%+03ddeg', frame, angle_deg);
            save(fullfile(lri_dir, [stem '.mat']), ...
                'env_img', 'bmode_img', 'frame', 'angle_idx', 'angle_deg', 'metadata');
            imwrite(png_img, fullfile(lri_dir, [stem '.png']));
        end
    end
end

function metadata = build_tracking_metadata(params, source)
    metadata.tracking_input = 'lri_env_frames';
    metadata.data_layout = '[z, x, angle, frame]';
    metadata.size_lri_env_frames = [params.Nz, params.Nx, params.n_angles, params.n_frames];
    metadata.is_compounded = false;
    metadata.source = source;
    metadata.run_id = get_optional_param(params, 'run_id', '');
    metadata.run_timestamp = get_optional_param(params, 'run_timestamp', '');
    metadata.output_root = get_optional_param(params, 'output_root', '');
    metadata.description = ...
        'Use lri_env_frames(:,:,angle_idx,frame_idx) for per-angle LRI speckle tracking.';

    metadata.x_grid_m = params.x_grid;
    metadata.z_grid_m = params.z_grid;
    metadata.dx_m = params.dx;
    metadata.dz_m = params.dz;
    metadata.angles_rad = params.angles;
    metadata.angles_deg = rad2deg(params.angles);
    metadata.PRF_Hz = params.PRF;
    metadata.dt_pulse_s = 1 / params.PRF;
    metadata.dt_frame_s = params.n_angles / params.PRF;
    metadata.array_order = 'row=z/depth, col=x/lateral';
    metadata.matlab_indexing = 'img = lri_env_frames(:,:,angle_idx,frame_idx)';
    metadata.not_tracking_inputs = {'hri_frames', 'env_frames', 'bmode_frames', 'compound_frames', 'videos'};
end

function value = get_optional_param(params, name, default_value)
    if isfield(params, name)
        value = params.(name);
    else
        value = default_value;
    end
end

function img = bmode_to_uint8(bmode, dynamic_range)
    scaled = 255 * (bmode + dynamic_range) / dynamic_range;
    scaled(bmode <= -dynamic_range) = 0;
    scaled(bmode >= 0) = 255;
    img = uint8(scaled);
end

function save_lri_montage(params, lri_bmode_frames, out_dir)
    figure('Name', 'Uncompounded LRI Angles', 'Position', [100, 540, 1200, 400]);
    for angle_idx = 1:params.n_angles
        subplot(1, params.n_angles, angle_idx);
        imagesc(params.x_grid * 1e3, params.z_grid * 1e3, lri_bmode_frames(:,:,angle_idx,1));
        colormap gray; colorbar;
        axis image;
        xlabel('Lateral [mm]');
        ylabel('Depth [mm]');
        title(sprintf('%.0f deg', rad2deg(params.angles(angle_idx))));
        clim([-params.dynamic_range 0]);
    end
    sgtitle('Uncompounded low-resolution images, frame 1');
    saveas(gcf, fullfile(out_dir, 'lri_angles_frame1.png'));
end
