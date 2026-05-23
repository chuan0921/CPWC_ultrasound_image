%% main_flow_simulation.m
% Fast image-domain CEUS flow phantom.
% Goal: reproduce the simulated B-mode appearance without Field II RF.

clearvars -except probe_name yscan_enabled yscan_positions_mm ...
    yscan_random_seeds yscan_dataset_name; close all; clc;

if ~exist('probe_name', 'var')
    probe_name = 'literature_l12_3v';
end

params = setup_parameters(probe_name);
if exist('yscan_enabled', 'var')
    params.yscan.enabled = logical(yscan_enabled);
end
if exist('yscan_positions_mm', 'var')
    params.yscan.y_positions_mm = yscan_positions_mm;
end
if exist('yscan_random_seeds', 'var')
    params.random_seeds = yscan_random_seeds;
    params.random_seed = params.random_seeds(1);
end
if exist('yscan_dataset_name', 'var')
    params.yscan.dataset_name = yscan_dataset_name;
end

proj_dir = fileparts(mfilename('fullpath'));
if isfield(params, 'yscan') && params.yscan.enabled
    run_yscan_batch(params, proj_dir);
    return;
end

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
    [lri_env_frames, lri_bmode_frames, lri_env_rows, lri_bmode_rows] = ...
        simulate_ceus_images(params);
    timing.total = toc(tic_seed);
    timing.seed = params.random_seed;
    fprintf('Image generation: %.2f s\n', timing.total);

    ground_truth = scatterer_medium('ground_truth', params);
    tracking_metadata = build_tracking_metadata(params, 'fast_image_domain');

    if isempty(lri_env_rows)
        save(fullfile(seed_dir, 'image_data.mat'), ...
            'lri_env_frames', 'lri_bmode_frames', ...
            'params', 'timing', 'tracking_metadata', '-v7.3');
    else
        save(fullfile(seed_dir, 'image_data.mat'), ...
            'lri_env_rows', 'lri_bmode_rows', ...
            'lri_env_frames', 'lri_bmode_frames', ...
            'params', 'timing', 'tracking_metadata', '-v7.3');
    end
    save(fullfile(seed_dir, 'ground_truth.mat'), 'ground_truth', '-v7.3');
    save(fullfile(seed_dir, 'tracking_metadata.mat'), 'tracking_metadata', '-v7.3');

    save_tracking_images(params, lri_env_frames, lri_bmode_frames, ...
        lri_env_rows, lri_bmode_rows, tracking_metadata, seed_dir);
    save_lri_montage(params, lri_bmode_frames, lri_bmode_rows, seed_dir);
    close all;

    seed_results{seed_idx}.seed = params.random_seed;
    seed_results{seed_idx}.out_dir = seed_dir;
    seed_results{seed_idx}.total_seconds = timing.total;
    fprintf('Data saved to: %s\n', seed_dir);
end

save(fullfile(run_out_dir, 'batch_summary.mat'), 'seed_results', 'params', '-v7.3');
fprintf('=== Done ===\n');

function run_yscan_batch(params, proj_dir)
    data_root = fullfile(proj_dir, params.yscan.output_root_name);
    if ~exist(data_root, 'dir'), mkdir(data_root); end

    run_out_dir = fullfile(data_root, params.yscan.dataset_name);
    if ~exist(run_out_dir, 'dir'), mkdir(run_out_dir); end

    params.run_id = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
    params.run_timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
    params.output_root = run_out_dir;

    fprintf('=== CPWC Y-scan Image Simulation Batch ===\n');
    fprintf('Probe: %s (%s)\n', params.probe.name, params.probe.label);
    fprintf('Dataset: %s\n', run_out_dir);
    fprintf('Y positions: %.2f to %.2f mm, spacing %.2f mm, count=%d\n', ...
        params.yscan.y_positions_mm(1), params.yscan.y_positions_mm(end), ...
        params.yscan.scan_spacing_mm, numel(params.yscan.y_positions_mm));
    fprintf('Grid: %d x %d pixels, frames/seed=%d, seeds=%d\n', ...
        params.Nz, params.Nx, params.n_frames, numel(params.random_seeds));

    seed_results = cell(numel(params.random_seeds), 1);
    for seed_idx = 1:numel(params.random_seeds)
        params.random_seed = params.random_seeds(seed_idx);
        rng(params.random_seed);

        seed_dir = fullfile(run_out_dir, sprintf('seed_%04d', params.random_seed));
        if ~exist(seed_dir, 'dir'), mkdir(seed_dir); end
        fprintf('--- Seed %d (%d/%d) ---\n', ...
            params.random_seed, seed_idx, numel(params.random_seeds));

        tic_seed = tic;
        yscan_data = simulate_yscan_images(params);
        timing.total = toc(tic_seed);
        timing.seed = params.random_seed;
        fprintf('Y-scan image generation: %.2f s\n', timing.total);

        yscan_metadata = build_yscan_metadata(params);
        save_json(fullfile(seed_dir, 'yscan_metadata.json'), yscan_metadata);
        save(fullfile(seed_dir, 'yscan_metadata.mat'), 'yscan_metadata', '-v7.3');

        n_rows = size(yscan_data, 1);
        for row_idx = 1:n_rows
            if n_rows > 1
                row_dir = fullfile(seed_dir, row_folder_name(row_idx));
            else
                row_dir = seed_dir;
            end
            if ~exist(row_dir, 'dir'), mkdir(row_dir); end

            template_dir = fullfile(row_dir, 'template_dictionary');
            if ~exist(template_dir, 'dir'), mkdir(template_dir); end

            for y_idx = 1:size(yscan_data, 2)
                y_data = yscan_data{row_idx, y_idx};
                y_dir = fullfile(row_dir, y_folder_name(y_data.true_y_mm));
                if ~exist(y_dir, 'dir'), mkdir(y_dir); end

                lri_env_frames = y_data.lri_env_frames;
                lri_bmode_frames = y_data.lri_bmode_frames;
                ground_truth = y_data.ground_truth;
                tracking_metadata = y_data.tracking_metadata;
                env_avg = y_data.env_avg;
                bmode_avg = y_data.bmode_avg;
                template_metadata = y_data.template_metadata;

                save(fullfile(y_dir, 'image_data.mat'), ...
                    'lri_env_frames', 'lri_bmode_frames', 'params', ...
                    'timing', 'tracking_metadata', '-v7.3');
                save(fullfile(y_dir, 'ground_truth.mat'), 'ground_truth', '-v7.3');
                save(fullfile(y_dir, 'tracking_metadata.mat'), ...
                    'tracking_metadata', '-v7.3');
                save(fullfile(y_dir, 'avg_image.mat'), ...
                    'env_avg', 'bmode_avg', 'template_metadata', '-v7.3');
                save_tracking_images(params, lri_env_frames, lri_bmode_frames, ...
                    [], [], tracking_metadata, y_dir);

                template_file = sprintf('%s_avg.mat', y_folder_name(y_data.true_y_mm));
                save(fullfile(template_dir, template_file), ...
                    'env_avg', 'bmode_avg', 'template_metadata', '-v7.3');
            end
        end

        close all;
        seed_results{seed_idx}.seed = params.random_seed;
        seed_results{seed_idx}.out_dir = seed_dir;
        seed_results{seed_idx}.total_seconds = timing.total;
        fprintf('Data saved to: %s\n', seed_dir);
    end

    save(fullfile(run_out_dir, 'batch_summary.mat'), ...
        'seed_results', 'params', '-v7.3');
    fprintf('=== Done ===\n');
end

function yscan_data = simulate_yscan_images(params)
    [X, Z] = meshgrid(params.x_grid, params.z_grid);
    flow_params = flow_domain_params(params);
    [X_flow, Z_flow] = meshgrid(flow_params.x_grid, flow_params.z_grid);
    [Z3, ~, Y3] = ndgrid(params.z_grid, params.x_grid, params.y_grid);
    [Z3_flow, ~, Y3_flow] = ndgrid(flow_params.z_grid, ...
        flow_params.x_grid, flow_params.y_grid);

    zc3 = yscan_center_z(params, Y3);
    zc3_flow = yscan_center_z(params, Y3_flow);
    r3 = sqrt(Y3.^2 + (Z3 - zc3).^2);
    r3_flow = sqrt(Y3_flow.^2 + (Z3_flow - zc3_flow).^2);
    lumen3 = r3 <= params.R;
    lumen3_flow = r3_flow <= params.R;
    wall3 = r3 > params.R & r3 <= params.R + params.wall_thickness;
    tissue3 = ~lumen3 & ~wall3;

    axial_psf = gaussian_kernel(params.lambda / 3, params.dz);
    lateral_psf = gaussian_kernel(params.lambda * 0.85, params.dx);
    elevation_psf = gaussian_kernel(params.lambda * 1.2, params.dy);

    tissue_iq = filtered_speckle_3d(params, axial_psf, lateral_psf, elevation_psf);
    flow_iq0 = filtered_speckle_3d(flow_params, ...
        axial_psf, lateral_psf, elevation_psf);
    wall_iq = filtered_speckle_3d(params, axial_psf, lateral_psf, elevation_psf);

    depth_gain = exp(-14 * (Z - params.z_min));
    lateral_gain = exp(-0.5 * (X / (0.42 * max(abs(params.x_grid)))).^2);
    sensitivity = depth_gain .* (0.55 + 0.45 * lateral_gain);

    tissue_iq = 0.28 * tissue_iq .* tissue3;
    flow_iq0 = 1.25 * flow_iq0 .* lumen3_flow;
    wall_iq = 2.6 * wall_iq .* wall3;

    vx_map3 = zeros(params.Nz, params.Nx, params.Ny);
    vx_map3(lumen3) = params.v0 * (1 - (r3(lumen3) / params.R).^2);

    y_positions_mm = params.yscan.y_positions_mm;
    n_ypos = numel(y_positions_mm);
    n_rows = yscan_row_count(params);
    yscan_data = cell(n_rows, n_ypos);
    for row_idx = 1:n_rows
        for y_idx = 1:n_ypos
            yscan_data{row_idx, y_idx}.true_y_mm = y_positions_mm(y_idx);
            yscan_data{row_idx, y_idx}.row_index = row_idx;
            yscan_data{row_idx, y_idx}.lri_env_frames = zeros( ...
                params.Nz, params.Nx, params.n_angles, params.n_frames);
        end
    end

    for frame = 1:params.n_frames
        for angle_idx = 1:params.n_angles
            pulse_idx = (frame - 1) * params.n_angles + (angle_idx - 1);
            displacement3 = vx_map3 * pulse_idx / params.PRF;
            flow_iq = shift_flow_volume_x(flow_iq0, displacement3, ...
                X, Z, X_flow, Z_flow);

            angle = params.angles(angle_idx);
            max_angle = max(abs(params.angles));
            if max_angle > 0
                angle_gain = 1 - 0.18 * abs(angle) / max_angle;
            else
                angle_gain = 1;
            end
            iq_volume = tissue_iq + flow_iq + wall_iq;

            for row_idx = 1:n_rows
                for y_idx = 1:n_ypos
                    y_m = y_positions_mm(y_idx) * 1e-3;
                    clean_iq = angle_gain * project_yscan_iq( ...
                        iq_volume, params, y_m, row_idx) ...
                        .* sensitivity .* yscan_shadow(params, y_m);
                    noisy_iq = clean_iq + complex_image_noise(size(X), clean_iq, params);
                    yscan_data{row_idx, y_idx}.lri_env_frames(:,:,angle_idx,frame) = ...
                        abs(noisy_iq);
                end
            end
        end
    end

    for row_idx = 1:n_rows
        for y_idx = 1:n_ypos
            y_m = y_positions_mm(y_idx) * 1e-3;
            y_data = yscan_data{row_idx, y_idx};
            reference = max(y_data.lri_env_frames(:));
            y_data.lri_bmode_frames = envelope_to_bmode( ...
                y_data.lri_env_frames, params.dynamic_range, reference);
            y_data.env_avg = mean(mean(y_data.lri_env_frames, 4), 3);
            y_data.bmode_avg = envelope_to_bmode( ...
                y_data.env_avg, params.dynamic_range, max(y_data.env_avg(:)));
            y_data.ground_truth = yscan_ground_truth(params, y_m, row_idx);
            y_data.tracking_metadata = build_yscan_tracking_metadata( ...
                params, y_m, row_idx);
            y_data.template_metadata = build_template_metadata(params, y_m, row_idx);
            yscan_data{row_idx, y_idx} = y_data;
        end
    end
end

function n_rows = yscan_row_count(params)
    if isfield(params.probe, 'n_y') && params.probe.n_y > 1
        n_rows = params.probe.n_y;
    else
        n_rows = 1;
    end
end

function img = project_yscan_iq(iq_volume, params, y_m, row_idx)
    if yscan_row_count(params) > 1 && ...
            isfield(params.probe, 'lr_geometry') && params.probe.lr_geometry
        weights = yscan_lr_slice_weights(params, y_m, row_idx);
        img = sum(iq_volume .* weights, 3);
        return;
    end

    weights = y_slice_weights(params, y_m);
    img = project_elevation(iq_volume, weights);
end

function weights = y_slice_weights(params, y_m)
    sigma = max(params.yscan.slice_sigma, params.dy);
    weights = exp(-0.5 * ((params.y_grid - y_m) / sigma).^2);
    weights = weights / sum(weights);
end

function weights = yscan_lr_slice_weights(params, y_m, row_idx)
    lr_distance = lr_distance_by_depth(params);
    center_sign = -1;
    if row_idx == 2
        center_sign = 1;
    end
    y_center = y_m + center_sign * lr_distance / 2;
    sigma = max(params.yscan.slice_sigma, params.dy);

    y = reshape(params.y_grid, 1, 1, []);
    y_center = reshape(y_center, [], 1, 1);
    weights = exp(-0.5 * ((y - y_center) / sigma).^2);
    weights = weights ./ max(sum(weights, 3), eps);
end

function shadow = yscan_shadow(params, y_m)
    [~, Z] = meshgrid(params.x_grid, params.z_grid);
    zc = yscan_center_z(params, y_m);
    r = sqrt(y_m.^2 + (Z - zc).^2);
    shadow = 1 - 0.18 * (r <= params.R);
end

function zc = yscan_center_z(params, y_m)
    z0 = params.yscan.z0_mm * 1e-3;
    slope = params.yscan.slope_mm_per_mm;
    zc = z0 + slope * y_m;
end

function ground_truth = yscan_ground_truth(params, y_m, row_idx)
    [X, Z] = meshgrid(params.x_grid, params.z_grid);
    zc = yscan_center_z(params, y_m);
    r = sqrt(y_m.^2 + (Z - zc).^2);
    inside = r <= params.R;

    vx_map = zeros(size(X));
    vx_map(inside) = params.v0 * (1 - (r(inside) / params.R).^2);

    dt_frame = params.n_angles / params.PRF;
    vessel = yscan_vessel_truth(params, y_m);

    ground_truth.vessel = vessel;
    ground_truth.row.index = row_idx;
    ground_truth.row.folder = row_folder_name(row_idx);
    ground_truth.row.effective_y_model = ...
        'y_eff(z)=scan_y+row_sign*LR(z)/2 for zipper LR geometry';
    ground_truth.flow.v0 = params.v0;
    ground_truth.flow.vx_map = vx_map;
    ground_truth.flow.vessel_mask = inside;
    ground_truth.flow.dt_frame = dt_frame;
    ground_truth.flow.dx_map = vx_map * dt_frame;
    ground_truth.flow.dx_max = max(vx_map(:)) * dt_frame;
    ground_truth.flow.flow_direction = 'x';
    ground_truth.grid.x = params.x_grid;
    ground_truth.grid.z = params.z_grid;
    ground_truth.grid.dx = params.dx;
    ground_truth.grid.dz = params.dz;
    ground_truth.grid.Nx = params.Nx;
    ground_truth.grid.Nz = params.Nz;
end

function vessel = yscan_vessel_truth(params, y_m)
    zc = yscan_center_z(params, y_m);
    R = params.R;
    R_outer = params.R + params.wall_thickness;
    y_abs = abs(y_m);

    vessel.true_y = y_m;
    vessel.true_y_mm = y_m * 1e3;
    vessel.center_z = zc;
    vessel.center_z_mm = zc * 1e3;
    vessel.radius = R;
    vessel.radius_mm = R * 1e3;
    vessel.wall_thickness = params.wall_thickness;
    vessel.wall_thickness_mm = params.wall_thickness * 1e3;

    if y_abs <= R
        lumen_half_z = sqrt(R^2 - y_abs^2);
        vessel.z_top = zc - lumen_half_z;
        vessel.z_bottom = zc + lumen_half_z;
        vessel.z_top_mm = vessel.z_top * 1e3;
        vessel.z_bottom_mm = vessel.z_bottom * 1e3;
    else
        vessel.z_top = NaN;
        vessel.z_bottom = NaN;
        vessel.z_top_mm = NaN;
        vessel.z_bottom_mm = NaN;
    end

    if y_abs <= R_outer
        outer_half_z = sqrt(R_outer^2 - y_abs^2);
        vessel.z_outer_top = zc - outer_half_z;
        vessel.z_outer_bottom = zc + outer_half_z;
        vessel.z_outer_top_mm = vessel.z_outer_top * 1e3;
        vessel.z_outer_bottom_mm = vessel.z_outer_bottom * 1e3;
    else
        vessel.z_outer_top = NaN;
        vessel.z_outer_bottom = NaN;
        vessel.z_outer_top_mm = NaN;
        vessel.z_outer_bottom_mm = NaN;
    end
end

function metadata = build_yscan_tracking_metadata(params, y_m, row_idx)
    metadata = build_tracking_metadata(params, 'fast_image_domain_yscan');
    metadata.tracking_input = 'lri_env_frames';
    metadata.data_layout = '[z, x, angle, frame]';
    metadata.matlab_indexing = 'img = lri_env_frames(:,:,angle_idx,frame_idx)';
    metadata.image_plane = 'xz';
    metadata.scan_axis = 'y';
    metadata.true_y_m = y_m;
    metadata.true_y_mm = y_m * 1e3;
    metadata.row_index = row_idx;
    metadata.row_folder = row_folder_name(row_idx);
    metadata.center_z_model = 'zc(y)=z0+slope*y';
    metadata.z0_mm = params.yscan.z0_mm;
    metadata.slope_mm_per_mm = params.yscan.slope_mm_per_mm;
    if yscan_row_count(params) > 1 && ...
            isfield(params.probe, 'lr_geometry') && params.probe.lr_geometry
        metadata.lr_depth_mm = params.probe.lr_depth_mm;
        metadata.lr_distance_mm = params.probe.lr_distance_mm;
        metadata.effective_y_model = ...
            'y_eff(z)=scan_y+row_sign*LR(z)/2';
    end
end

function metadata = build_template_metadata(params, y_m, row_idx)
    vessel = yscan_vessel_truth(params, y_m);
    metadata.true_y_mm = vessel.true_y_mm;
    metadata.row_index = row_idx;
    metadata.row_folder = row_folder_name(row_idx);
    metadata.z_top_truth_mm = vessel.z_top_mm;
    metadata.z_bottom_truth_mm = vessel.z_bottom_mm;
    metadata.z_outer_top_truth_mm = vessel.z_outer_top_mm;
    metadata.z_outer_bottom_truth_mm = vessel.z_outer_bottom_mm;
    metadata.center_z_truth_mm = vessel.center_z_mm;
    metadata.radius_truth_mm = vessel.radius_mm;
    metadata.wall_thickness_mm = vessel.wall_thickness_mm;
    metadata.center_z_model = 'zc(y)=z0+slope*y';
    metadata.z0_mm = params.yscan.z0_mm;
    metadata.slope_mm_per_mm = params.yscan.slope_mm_per_mm;
    if yscan_row_count(params) > 1 && ...
            isfield(params.probe, 'lr_geometry') && params.probe.lr_geometry
        metadata.lr_depth_mm = params.probe.lr_depth_mm;
        metadata.lr_distance_mm = params.probe.lr_distance_mm;
        metadata.effective_y_model = ...
            'y_eff(z)=scan_y+row_sign*LR(z)/2';
    end
end

function metadata = build_yscan_metadata(params)
    metadata.scan_axis = 'y';
    metadata.y_positions_mm = params.yscan.y_positions_mm;
    metadata.scan_spacing_mm = params.yscan.scan_spacing_mm;
    metadata.vessel_radius_mm = params.R * 1e3;
    metadata.vessel_inner_diameter_mm = 2 * params.R * 1e3;
    metadata.wall_thickness_mm = params.wall_thickness * 1e3;
    metadata.center_z_model = 'zc(y)=z0+slope*y';
    metadata.z0_mm = params.yscan.z0_mm;
    metadata.slope_mm_per_mm = params.yscan.slope_mm_per_mm;
    metadata.flow_direction = 'x';
    metadata.image_plane = 'xz';
    metadata.probe_name = params.probe.name;
    metadata.f0_MHz = params.f0 / 1e6;
    metadata.PRF_Hz = params.PRF;
    metadata.n_angles = params.n_angles;
    metadata.n_frames = params.n_frames;
    metadata.dx_mm = params.dx * 1e3;
    metadata.dz_mm = params.dz * 1e3;
    metadata.slice_sigma_mm = params.yscan.slice_sigma * 1e3;
    metadata.row_count = yscan_row_count(params);
    if metadata.row_count > 1
        metadata.row_folders = arrayfun(@row_folder_name, ...
            1:metadata.row_count, 'UniformOutput', false);
        metadata.folder_layout = 'seed/rowN/y_*.mm/';
    else
        metadata.folder_layout = 'seed/y_*.mm/';
    end
end

function name = y_folder_name(y_mm)
    if y_mm < 0
        sign_char = '-';
    else
        sign_char = '+';
    end
    name = sprintf('y_%s%0.2fmm', sign_char, abs(y_mm));
end

function name = row_folder_name(row_idx)
    name = sprintf('row%d', row_idx);
end

function save_json(filename, data)
    text = jsonencode(data, PrettyPrint=true);
    fid = fopen(filename, 'w');
    if fid < 0
        error('Unable to open JSON output: %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, '%s\n', text);
end

function [lri_env_frames, lri_bmode_frames, lri_env_rows, lri_bmode_rows] = ...
        simulate_ceus_images(params)
    [X, Z] = meshgrid(params.x_grid, params.z_grid);
    flow_params = flow_domain_params(params);
    [X_flow, Z_flow] = meshgrid(flow_params.x_grid, flow_params.z_grid);
    [Z3, ~, Y3] = ndgrid(params.z_grid, params.x_grid, params.y_grid);
    [Z3_flow, ~, Y3_flow] = ndgrid(flow_params.z_grid, ...
        flow_params.x_grid, flow_params.y_grid);

    z0 = params.vessel_center_z;
    r3 = sqrt(Y3.^2 + (Z3 - z0).^2);
    r3_flow = sqrt(Y3_flow.^2 + (Z3_flow - z0).^2);
    lumen3 = r3 <= params.R;
    lumen3_flow = r3_flow <= params.R;
    wall3 = r3 > params.R & r3 <= params.R + params.wall_thickness;
    tissue3 = ~lumen3 & ~wall3;

    central_y_idx = nearest_index(params.y_grid, 0);
    lumen = lumen3(:,:,central_y_idx);

    % Image-domain point spread function. This approximates compounded plane
    % wave resolution without simulating channel RF.
    axial_psf = gaussian_kernel(params.lambda / 3, params.dz);
    lateral_psf = gaussian_kernel(params.lambda * 0.85, params.dx);
    elevation_psf = gaussian_kernel(params.lambda * 1.2, params.dy);

    tissue_iq = filtered_speckle_3d(params, axial_psf, lateral_psf, elevation_psf);
    flow_iq0 = filtered_speckle_3d(flow_params, ...
        axial_psf, lateral_psf, elevation_psf);
    wall_iq = filtered_speckle_3d(params, axial_psf, lateral_psf, elevation_psf);

    has_rows = isfield(params.probe, 'n_y') && params.probe.n_y > 1;
    if has_rows
        n_rows = params.probe.n_y;
    else
        n_rows = 0;
    end

    depth_gain = exp(-14 * (Z - params.z_min));
    lateral_gain = exp(-0.5 * (X / (0.42 * max(abs(params.x_grid)))).^2);
    sensitivity = depth_gain .* (0.55 + 0.45 * lateral_gain);

    tissue_iq = 0.28 * tissue_iq .* tissue3;
    flow_iq0 = 1.25 * flow_iq0 .* lumen3_flow;
    wall_iq = 2.6 * wall_iq .* wall3;

    vx_map3 = zeros(params.Nz, params.Nx, params.Ny);
    vx_map3(lumen3) = params.v0 * (1 - (r3(lumen3) / params.R).^2);

    lri_env_frames = zeros(params.Nz, params.Nx, params.n_angles, params.n_frames);
    lri_bmode_frames = zeros(params.Nz, params.Nx, params.n_angles, params.n_frames);
    if has_rows
        lri_env_rows = zeros(params.Nz, params.Nx, n_rows, params.n_angles, params.n_frames);
        lri_bmode_rows = zeros(params.Nz, params.Nx, n_rows, params.n_angles, params.n_frames);
    else
        lri_env_rows = [];
        lri_bmode_rows = [];
    end

    for frame = 1:params.n_frames
        for angle_idx = 1:params.n_angles
            pulse_idx = (frame - 1) * params.n_angles + (angle_idx - 1);
            displacement3 = vx_map3 * pulse_idx / params.PRF;
            flow_iq = shift_flow_volume_x(flow_iq0, displacement3, ...
                X, Z, X_flow, Z_flow);

            angle = params.angles(angle_idx);
            max_angle = max(abs(params.angles));
            if max_angle > 0
                angle_gain = 1 - 0.18 * abs(angle) / max_angle;
            else
                angle_gain = 1;
            end
            shadow = 1 - 0.18 * lumen;
            iq_volume = tissue_iq + flow_iq + wall_iq;

            if has_rows
                row_iq_stack = zeros(params.Nz, params.Nx, n_rows);
                for row = 1:n_rows
                    row_gain = 1 + 0.03 * (row - (n_rows + 1) / 2);
                    row_iq = project_row_iq(iq_volume, params, row);
                    row_iq = row_gain * angle_gain * row_iq .* sensitivity .* shadow;
                    row_iq_stack(:,:,row) = row_iq;
                end
                clean_lri_iq = mean(row_iq_stack, 3);
                lri_iq = clean_lri_iq + complex_image_noise(size(X), clean_lri_iq, params);
            else
                aperture_weights = aperture_elevation_weights(params);
                clean_lri_iq = angle_gain * project_elevation(iq_volume, aperture_weights) ...
                    .* sensitivity .* shadow;
                lri_iq = clean_lri_iq + complex_image_noise(size(X), clean_lri_iq, params);
            end

            lri_env_frames(:,:,angle_idx,frame) = abs(lri_iq);

            for row = 1:n_rows
                row_iq = row_iq_stack(:,:,row) + ...
                    complex_image_noise(size(X), row_iq_stack(:,:,row), params);

                lri_env_rows(:,:,row,angle_idx,frame) = abs(row_iq);
            end
        end
    end

    bmode_reference = bmode_reference_value(lri_env_frames, lri_env_rows, params);
    lri_bmode_frames = envelope_to_bmode(lri_env_frames, ...
        params.dynamic_range, bmode_reference);
    if has_rows
        lri_bmode_rows = envelope_to_bmode(lri_env_rows, ...
            params.dynamic_range, bmode_reference);
    end

    expected_size = [params.Nz, params.Nx, params.n_angles, params.n_frames];
    assert(isequal(size(lri_env_frames), expected_size), ...
        'lri_env_frames must be [z, x, angle, frame].');
    if has_rows
        expected_row_size = [params.Nz, params.Nx, n_rows, ...
            params.n_angles, params.n_frames];
        assert(isequal(size(lri_env_rows), expected_row_size), ...
            'lri_env_rows must be [z, x, row, angle, frame].');
    end
end

function idx = nearest_index(grid, value)
    [~, idx] = min(abs(grid - value));
end

function flow_params = flow_domain_params(params)
    flow_params = params;
    mode = get_optional_param(params, 'flow_domain_mode', 'fov');
    switch lower(mode)
        case 'extended_crop'
            flow_params.x_grid = params.flow_x_grid;
            flow_params.Nx = params.Nx_flow;
        case 'fov'
        otherwise
            error('Unknown flow_domain_mode: %s', mode);
    end
end

function iq = filtered_speckle_3d(params, axial_psf, lateral_psf, elevation_psf)
    iq = complex(randn(params.Nz, params.Nx, params.Ny), ...
        randn(params.Nz, params.Nx, params.Ny));
    iq = convn(iq, reshape(axial_psf, [], 1, 1), 'same');
    iq = convn(iq, reshape(lateral_psf, 1, [], 1), 'same');
    iq = convn(iq, reshape(elevation_psf, 1, 1, []), 'same');
end

function flow_iq = shift_flow_volume_x(flow_iq0, displacement3, X, Z, X_flow, Z_flow)
    flow_iq = complex(zeros(size(displacement3)));
    for y_idx = 1:size(flow_iq0, 3)
        flow_iq(:,:,y_idx) = interp2(X_flow, Z_flow, flow_iq0(:,:,y_idx), ...
            X - displacement3(:,:,y_idx), Z, 'linear', 0);
    end
end

function img = project_elevation(iq_volume, weights)
    weights = reshape(weights, 1, 1, []);
    img = sum(iq_volume .* weights, 3);
end

function img = project_row_iq(iq_volume, params, row)
    if isfield(params.probe, 'lr_geometry') && params.probe.lr_geometry
        img = project_lr_geometry(iq_volume, params, row);
        return;
    end

    row_weights = row_elevation_weights(params, row);
    img = project_elevation(iq_volume, row_weights);
end

function img = project_lr_geometry(iq_volume, params, row)
    lr_distance = lr_distance_by_depth(params);
    center_sign = -1;
    if row == 2
        center_sign = 1;
    end
    y_center = center_sign * lr_distance / 2;
    sigma = max(params.probe.lr_sigma, params.dy);

    y = reshape(params.y_grid, 1, 1, []);
    y_center = reshape(y_center, [], 1, 1);
    weights = exp(-0.5 * ((y - y_center) / sigma).^2);
    weights = weights ./ max(sum(weights, 3), eps);
    img = sum(iq_volume .* weights, 3);
end

function lr_distance = lr_distance_by_depth(params)
    depth_m = params.probe.lr_depth_mm * 1e-3;
    distance_m = params.probe.lr_distance_mm * 1e-3;
    lr_distance = interp1(depth_m, distance_m, params.z_grid, ...
        'linear', 'extrap')';
end

function weights = aperture_elevation_weights(params)
    if isfield(params.probe, 'element_data')
        y_vertices = params.probe.element_data(:, [3 6 9 12]);
        y_min = min(y_vertices(:));
        y_max = max(y_vertices(:));
        weights = soft_rect_weights(params.y_grid, y_min, y_max, params.dy);
    else
        sigma = max(params.elevation_beam_sigma, params.lambda);
        weights = exp(-0.5 * (params.y_grid / sigma).^2);
    end
    weights = weights / sum(weights);
end

function noise = complex_image_noise(image_size, clean_iq, params)
    mode = get_optional_param(params, 'noise_mode', 'floor');
    switch lower(mode)
        case 'snr'
            signal_rms = sqrt(mean(abs(clean_iq(:)).^2));
            complex_noise_rms = signal_rms / 10^(params.SNR_dB / 20);
            component_sigma = complex_noise_rms / sqrt(2);
            noise = component_sigma * complex(randn(image_size), randn(image_size));
        case 'floor'
            noise = params.image_noise_floor * complex(randn(image_size), randn(image_size));
        otherwise
            error('Unknown noise_mode: %s', mode);
    end
end

function weights = row_elevation_weights(params, row)
    element_rows = row:params.probe.n_y:params.probe.n_elements;
    y_vertices = params.probe.element_data(element_rows, [3 6 9 12]);
    y_min = min(y_vertices(:));
    y_max = max(y_vertices(:));
    weights = soft_rect_weights(params.y_grid, y_min, y_max, params.dy);
    weights = weights / sum(weights);
end

function weights = soft_rect_weights(y_grid, y_min, y_max, edge_width)
    left = 1 ./ (1 + exp(-(y_grid - y_min) / edge_width));
    right = 1 ./ (1 + exp((y_grid - y_max) / edge_width));
    weights = left .* right;
    if sum(weights) == 0
        weights(:) = 1;
    end
end

function kernel = gaussian_kernel(sigma_m, spacing_m)
    radius = max(1, ceil(4 * sigma_m / spacing_m));
    x = (-radius:radius) * spacing_m;
    kernel = exp(-0.5 * (x / sigma_m).^2);
    kernel = kernel / sum(kernel);
end

function reference = bmode_reference_value(lri_env_frames, lri_env_rows, params)
    mode = get_optional_param(params, 'bmode_reference_mode', 'seed_max');
    switch lower(mode)
        case 'seed_max'
            if isempty(lri_env_rows)
                reference = max(lri_env_frames(:));
            else
                reference = max([lri_env_frames(:); lri_env_rows(:)]);
            end
        otherwise
            error('Unknown bmode_reference_mode: %s', mode);
    end
end

function bmode = envelope_to_bmode(envelope, dynamic_range, reference)
    bmode = 20 * log10(max(envelope, eps) / max(reference, eps));
    bmode = max(bmode, -dynamic_range);
end

function save_tracking_images(params, lri_env_frames, lri_bmode_frames, ...
        lri_env_rows, lri_bmode_rows, metadata, out_dir)
    lri_dir = fullfile(out_dir, 'lri_frames');
    if ~exist(lri_dir, 'dir'), mkdir(lri_dir); end

    if isempty(lri_env_rows)
        for frame = 1:params.n_frames
            for angle_idx = 1:params.n_angles
                angle_deg = round(rad2deg(params.angles(angle_idx)));
                env_img = lri_env_frames(:,:,angle_idx,frame);
                bmode_img = lri_bmode_frames(:,:,angle_idx,frame);
                png_img = bmode_to_uint8(bmode_img, params.dynamic_range);

                stem = sprintf('frame_%03d_angle_%+03ddeg', frame, angle_deg);
                save(fullfile(lri_dir, [stem '.mat']), ...
                    'env_img', 'bmode_img', 'frame', 'angle_idx', ...
                    'angle_deg', 'metadata');
                imwrite(png_img, fullfile(lri_dir, [stem '.png']));
            end
        end
        return;
    end

    for row = 1:size(lri_env_rows, 3)
        row_dir = fullfile(lri_dir, sprintf('row_%02d', row));
        if ~exist(row_dir, 'dir'), mkdir(row_dir); end

        for frame = 1:params.n_frames
            for angle_idx = 1:params.n_angles
                angle_deg = round(rad2deg(params.angles(angle_idx)));
                env_img = lri_env_rows(:,:,row,angle_idx,frame);
                bmode_img = lri_bmode_rows(:,:,row,angle_idx,frame);
                png_img = bmode_to_uint8(bmode_img, params.dynamic_range);

                stem = sprintf('frame_%03d_row_%02d_angle_%+03ddeg', ...
                    frame, row, angle_deg);
                save(fullfile(row_dir, [stem '.mat']), ...
                    'env_img', 'bmode_img', 'frame', 'row', 'angle_idx', ...
                    'angle_deg', 'metadata');
                imwrite(png_img, fullfile(row_dir, [stem '.png']));
            end
        end
    end
end

function metadata = build_tracking_metadata(params, source)
    has_rows = isfield(params.probe, 'n_y') && params.probe.n_y > 1;
    if has_rows
        metadata.tracking_input = 'lri_env_rows';
        metadata.data_layout = '[z, x, row, angle, frame]';
        metadata.size_lri_env_rows = [params.Nz, params.Nx, ...
            params.probe.n_y, params.n_angles, params.n_frames];
        metadata.row_meaning = 'row 1 = odd elements, row 2 = even elements';
    else
        metadata.tracking_input = 'lri_env_frames';
        metadata.data_layout = '[z, x, angle, frame]';
    end
    metadata.size_lri_env_frames = [params.Nz, params.Nx, params.n_angles, params.n_frames];
    metadata.is_compounded = false;
    metadata.source = source;
    metadata.run_id = get_optional_param(params, 'run_id', '');
    metadata.run_timestamp = get_optional_param(params, 'run_timestamp', '');
    metadata.output_root = get_optional_param(params, 'output_root', '');
    if has_rows
        metadata.description = ...
            'Use lri_env_rows(:,:,row_idx,angle_idx,frame_idx) for zipper row-specific LRI tracking.';
    else
        metadata.description = ...
            'Use lri_env_frames(:,:,angle_idx,frame_idx) for per-angle LRI speckle tracking.';
    end

    metadata.x_grid_m = params.x_grid;
    metadata.z_grid_m = params.z_grid;
    metadata.y_grid_m = params.y_grid;
    metadata.dx_m = params.dx;
    metadata.dy_m = params.dy;
    metadata.dz_m = params.dz;
    metadata.angles_rad = params.angles;
    metadata.angles_deg = rad2deg(params.angles);
    metadata.PRF_Hz = params.PRF;
    metadata.dt_pulse_s = 1 / params.PRF;
    metadata.dt_frame_s = params.n_angles / params.PRF;
    metadata.array_order = 'row=z/depth, col=x/lateral';
    if has_rows
        metadata.matlab_indexing = ...
            'img = lri_env_rows(:,:,row_idx,angle_idx,frame_idx)';
        metadata.row_y_centers_m = row_y_centers(params);
        metadata.row_projection = ...
            'Each row image is an elevation-weighted x-z projection from the same 3D y-z vessel phantom.';
        if isfield(params.probe, 'lr_geometry') && params.probe.lr_geometry
            metadata.row_projection = ...
                'Rows use depth-dependent zipper LR geometry in y/elevation.';
            metadata.lr_depth_mm = params.probe.lr_depth_mm;
            metadata.lr_distance_mm = params.probe.lr_distance_mm;
            metadata.lr_sigma_m = params.probe.lr_sigma;
        end
    else
        metadata.matlab_indexing = 'img = lri_env_frames(:,:,angle_idx,frame_idx)';
    end
    metadata.not_tracking_inputs = {'hri_frames', 'env_frames', 'bmode_frames', 'compound_frames', 'videos'};
end

function centers = row_y_centers(params)
    centers = zeros(1, params.probe.n_y);
    for row = 1:params.probe.n_y
        element_rows = row:params.probe.n_y:params.probe.n_elements;
        y_vertices = params.probe.element_data(element_rows, [3 6 9 12]);
        centers(row) = mean([min(y_vertices(:)), max(y_vertices(:))]);
    end
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

function save_lri_montage(params, lri_bmode_frames, lri_bmode_rows, out_dir)
    if isempty(lri_bmode_rows)
        n_rows = 1;
    else
        n_rows = size(lri_bmode_rows, 3);
    end

    figure('Name', 'Uncompounded LRI Angles', ...
        'Position', [100, 540, 1200, 260 * n_rows]);
    for row = 1:n_rows
        for angle_idx = 1:params.n_angles
            subplot(n_rows, params.n_angles, (row - 1) * params.n_angles + angle_idx);
            if isempty(lri_bmode_rows)
                img = lri_bmode_frames(:,:,angle_idx,1);
                row_title = '';
            else
                img = lri_bmode_rows(:,:,row,angle_idx,1);
                row_title = sprintf('row %d, ', row);
            end
            imagesc(params.x_grid * 1e3, params.z_grid * 1e3, img);
            colormap gray; colorbar;
            axis image;
            xlabel('Lateral [mm]');
            ylabel('Depth [mm]');
            title(sprintf('%s%.0f deg', row_title, rad2deg(params.angles(angle_idx))));
            clim([-params.dynamic_range 0]);
        end
    end
    sgtitle('Uncompounded low-resolution images, frame 1');
    saveas(gcf, fullfile(out_dir, 'lri_angles_frame1.png'));
end
