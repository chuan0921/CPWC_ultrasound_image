%% main_flow_simulation.m
% Fast image-domain CEUS flow phantom.
% Goal: reproduce the simulated B-mode appearance without Field II RF.

clearvars -except probe_name; close all; clc;

if ~exist('probe_name', 'var')
    probe_name = 'literature_l12_3v';
end

params = setup_parameters(probe_name);

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

function [lri_env_frames, lri_bmode_frames, lri_env_rows, lri_bmode_rows] = ...
        simulate_ceus_images(params)
    [X, Z] = meshgrid(params.x_grid, params.z_grid);
    [Z3, ~, Y3] = ndgrid(params.z_grid, params.x_grid, params.y_grid);

    z0 = params.vessel_center_z;
    r3 = sqrt(Y3.^2 + (Z3 - z0).^2);
    lumen3 = r3 <= params.R;
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
    flow_iq0 = filtered_speckle_3d(params, axial_psf, lateral_psf, elevation_psf);
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
    flow_iq0 = 1.25 * flow_iq0 .* lumen3;
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
            flow_iq = shift_flow_volume_x(flow_iq0, displacement3, X, Z);

            angle = params.angles(angle_idx);
            max_angle = max(abs(params.angles));
            if max_angle > 0
                angle_gain = 1 - 0.18 * abs(angle) / max_angle;
            else
                angle_gain = 1;
            end
            angle_noise = params.image_noise_floor * complex(randn(size(X)), randn(size(X)));

            shadow = 1 - 0.18 * lumen;
            iq_volume = tissue_iq + flow_iq + wall_iq;

            if has_rows
                row_iq_stack = zeros(params.Nz, params.Nx, n_rows);
                for row = 1:n_rows
                    row_weights = row_elevation_weights(params, row);
                    row_gain = 1 + 0.03 * (row - (n_rows + 1) / 2);
                    row_iq = project_elevation(iq_volume, row_weights);
                    row_iq = row_gain * angle_gain * row_iq .* sensitivity .* shadow;
                    row_iq_stack(:,:,row) = row_iq;
                end
                lri_iq = mean(row_iq_stack, 3) + angle_noise;
            else
                aperture_weights = aperture_elevation_weights(params);
                lri_iq = angle_gain * project_elevation(iq_volume, aperture_weights) ...
                    .* sensitivity .* shadow;
                lri_iq = lri_iq + angle_noise;
            end

            lri_env_frames(:,:,angle_idx,frame) = abs(lri_iq);
            lri_bmode_frames(:,:,angle_idx,frame) = ...
                envelope_to_bmode(lri_env_frames(:,:,angle_idx,frame), params.dynamic_range);

            for row = 1:n_rows
                row_noise = params.image_noise_floor * ...
                    complex(randn(size(X)), randn(size(X)));
                row_iq = row_iq_stack(:,:,row) + row_noise;

                lri_env_rows(:,:,row,angle_idx,frame) = abs(row_iq);
                lri_bmode_rows(:,:,row,angle_idx,frame) = envelope_to_bmode( ...
                    lri_env_rows(:,:,row,angle_idx,frame), params.dynamic_range);
            end
        end
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

function iq = filtered_speckle_3d(params, axial_psf, lateral_psf, elevation_psf)
    iq = complex(randn(params.Nz, params.Nx, params.Ny), ...
        randn(params.Nz, params.Nx, params.Ny));
    iq = convn(iq, reshape(axial_psf, [], 1, 1), 'same');
    iq = convn(iq, reshape(lateral_psf, 1, [], 1), 'same');
    iq = convn(iq, reshape(elevation_psf, 1, 1, []), 'same');
end

function flow_iq = shift_flow_volume_x(flow_iq0, displacement3, X, Z)
    flow_iq = complex(zeros(size(flow_iq0)));
    for y_idx = 1:size(flow_iq0, 3)
        flow_iq(:,:,y_idx) = interp2(X, Z, flow_iq0(:,:,y_idx), ...
            X - displacement3(:,:,y_idx), Z, 'linear', 0);
    end
end

function img = project_elevation(iq_volume, weights)
    weights = reshape(weights, 1, 1, []);
    img = sum(iq_volume .* weights, 3);
end

function weights = aperture_elevation_weights(params)
    if isfield(params.probe, 'element_data')
        y_vertices = params.probe.element_data(:, [3 6 9 12]);
        y_min = min(y_vertices(:));
        y_max = max(y_vertices(:));
        weights = soft_rect_weights(params.y_grid, y_min, y_max, params.dy);
    else
        sigma = max(params.R / 2, params.lambda);
        weights = exp(-0.5 * (params.y_grid / sigma).^2);
    end
    weights = weights / sum(weights);
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

function bmode = envelope_to_bmode(envelope, dynamic_range)
    env_max = max(envelope(:));
    bmode = 20 * log10(max(envelope, eps) / max(env_max, eps));
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
