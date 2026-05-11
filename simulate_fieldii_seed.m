function result = simulate_fieldii_seed(params, out_dir, fieldii_dir)
% simulate_fieldii_seed - Run one seed; angles parallelized via parfor.

    addpath(fieldii_dir);
    rng(params.random_seed);
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    fprintf('[seed %d] start\n', params.random_seed);
    tic_seed = tic;

    [flow_pos0, flow_amp, static_pos, static_amp] = scatterer_medium('generate', params);

    R = params.R;
    z_c = params.vessel_center_z;
    r0 = sqrt(flow_pos0(:,2).^2 + (flow_pos0(:,3) - z_c).^2);
    vx = params.v0 * max(0, 1 - (r0 / R).^2);

    dt = 1 / params.PRF;
    n_angles = params.n_angles;
    n_frames = params.n_frames;
    angles_array = params.angles;
    seed_id = params.random_seed;

    angle_results = cell(n_angles, 1);
    parfor a = 1:n_angles
        angle_results{a} = simulate_one_angle( ...
            params, fieldii_dir, ...
            static_pos, static_amp, ...
            flow_pos0, flow_amp, vx, ...
            angles_array(a), a, dt, n_frames, seed_id);
    end

    lri_frames = zeros(params.Nz, params.Nx, n_angles, n_frames);
    for a = 1:n_angles
        lri_frames(:,:,a,:) = angle_results{a}.lri;
    end

    lri_env_frames = zeros(params.Nz, params.Nx, n_angles, n_frames);
    lri_bmode_frames = zeros(params.Nz, params.Nx, n_angles, n_frames);
    for frame = 1:n_frames
        for a = 1:n_angles
            lri_rf = lri_frames(:,:,a,frame);
            lri_env_frames(:,:,a,frame) = abs(hilbert(lri_rf));
            lri_bmode_frames(:,:,a,frame) = envelope_to_bmode( ...
                lri_env_frames(:,:,a,frame), params.dynamic_range);
        end
    end
    expected_size = [params.Nz, params.Nx, params.n_angles, params.n_frames];
    assert(isequal(size(lri_env_frames), expected_size), ...
        'lri_env_frames must be [z, x, angle, frame].');

    timing.static_rf = zeros(n_angles, 1);
    timing.flow_rf = zeros(n_frames, n_angles);
    timing.das = zeros(n_frames, n_angles);
    for a = 1:n_angles
        timing.static_rf(a) = angle_results{a}.timing.static;
        timing.flow_rf(:, a) = angle_results{a}.timing.flow;
        timing.das(:, a) = angle_results{a}.timing.das;
    end
    timing.total = toc(tic_seed);

    ground_truth = scatterer_medium('ground_truth', params);
    tracking_metadata = build_tracking_metadata(params, 'fieldii_das');

    save(fullfile(out_dir, 'image_data.mat'), ...
        'lri_env_frames', 'lri_bmode_frames', ...
        'params', 'timing', 'tracking_metadata', '-v7.3');
    save(fullfile(out_dir, 'ground_truth.mat'), 'ground_truth', '-v7.3');
    save(fullfile(out_dir, 'tracking_metadata.mat'), 'tracking_metadata', '-v7.3');
    save_lri_pngs(params, lri_bmode_frames, out_dir);

    result.seed = seed_id;
    result.out_dir = out_dir;
    result.total_seconds = timing.total;
    result.static_rf_seconds = sum(timing.static_rf(:));
    result.flow_rf_seconds = sum(timing.flow_rf(:));
    result.das_seconds = sum(timing.das(:));

    fprintf('[seed %d] done: %.2f s\n', seed_id, timing.total);
end

function out = simulate_one_angle(params, fieldii_dir, ...
        static_pos, static_amp, flow_pos0, flow_amp, vx, ...
        angle, angle_idx, dt, n_frames, seed_id)
% Worker entry: runs 1 angle's full pipeline (static + n_frames flow + DAS).

    addpath(fieldii_dir);

    field_init(0);
    cleanup_obj = onCleanup(@() cleanup_fieldii());
    set_field('fs', params.fs);
    set_field('c', params.c);

    [Tx, Rx] = create_probe(params);
    probe_cleanup = onCleanup(@() cleanup_probe(Tx, Rx));

    noise_stream = RandStream('twister', 'Seed', seed_id * 1000 + angle_idx);

    tic_s = tic;
    [static_rf, static_t0] = simulate_rf_planewave_clean( ...
        Tx, Rx, static_pos, static_amp, angle, params);
    static_time = toc(tic_s);
    fprintf('[seed %d angle %d] static: %.2f s\n', seed_id, angle_idx, static_time);

    n_angles = params.n_angles;
    half_len = params.vessel_length / 2;
    lri_stack = zeros(params.Nz, params.Nx, n_frames);
    flow_times = zeros(n_frames, 1);
    das_times = zeros(n_frames, 1);

    for frame = 1:n_frames
        t = ((frame - 1) * n_angles + angle_idx) * dt;
        x_t = mod(flow_pos0(:,1) + vx * t + half_len, 2*half_len) - half_len;
        flow_pos_t = [x_t, flow_pos0(:,2), flow_pos0(:,3)];

        tic_f = tic;
        [flow_rf, flow_t0] = simulate_rf_planewave_clean( ...
            Tx, Rx, flow_pos_t, flow_amp, angle, params);
        [rf_data, t_start] = add_rf_signals( ...
            static_rf, static_t0, flow_rf, flow_t0, params.fs);
        rf_data = add_rf_noise_stream(rf_data, params.SNR_dB, noise_stream);
        flow_times(frame) = toc(tic_f);

        tic_d = tic;
        lri_stack(:,:,frame) = das_beamform_planewave( ...
            rf_data, t_start, angle, params);
        das_times(frame) = toc(tic_d);

        fprintf('[seed %d angle %d] frame %d/%d: flow %.2f s, DAS %.2f s\n', ...
            seed_id, angle_idx, frame, n_frames, flow_times(frame), das_times(frame));
    end

    out.lri = lri_stack;
    out.timing.static = static_time;
    out.timing.flow = flow_times;
    out.timing.das = das_times;
end

function cleanup_probe(Tx, Rx)
    try
        xdc_free(Tx);
    catch
    end
    try
        xdc_free(Rx);
    catch
    end
end

function cleanup_fieldii()
    try
        field_end;
    catch
    end
end

function [rf_data, t_start] = simulate_rf_planewave_clean(Tx, Rx, positions, amplitudes, angle, params)
    elem_x = ((0:params.n_elements-1) - (params.n_elements-1)/2) * params.pitch;
    delays = -elem_x * sin(angle) / params.c;
    delays = delays - min(delays);

    xdc_focus_times(Tx, 0, delays);
    xdc_focus_times(Rx, 0, zeros(1, params.n_elements));

    [rf_data, t_start] = calc_scat_multi(Tx, Rx, positions, amplitudes);
end

function [rf_sum, t_sum] = add_rf_signals(rf_a, t_a, rf_b, t_b, fs)
    if isempty(rf_a)
        rf_sum = rf_b;
        t_sum = t_b;
        return;
    end
    if isempty(rf_b)
        rf_sum = rf_a;
        t_sum = t_a;
        return;
    end

    n_elem = max(size(rf_a, 2), size(rf_b, 2));
    t_sum = min(t_a, t_b);

    off_a = round((t_a - t_sum) * fs);
    off_b = round((t_b - t_sum) * fs);
    n_rows = max(off_a + size(rf_a, 1), off_b + size(rf_b, 1));

    rf_sum = zeros(n_rows, n_elem);
    rows_a = off_a + (1:size(rf_a, 1));
    rows_b = off_b + (1:size(rf_b, 1));

    rf_sum(rows_a, 1:size(rf_a, 2)) = rf_sum(rows_a, 1:size(rf_a, 2)) + rf_a;
    rf_sum(rows_b, 1:size(rf_b, 2)) = rf_sum(rows_b, 1:size(rf_b, 2)) + rf_b;
end

function rf_data = add_rf_noise_stream(rf_data, snr_db, stream)
    sig_power = mean(rf_data(:).^2);
    if sig_power > 0
        noise_power = sig_power / (10^(snr_db / 10));
        rf_data = rf_data + sqrt(noise_power) * randn(stream, size(rf_data));
    end
end

function [Tx, Rx] = create_probe(params)
    if isfield(params, 'probe') && ~strcmp(params.probe.geometry, 'linear_array')
        error('create_probe only supports linear_array Field II xdc. Selected probe is %s.', ...
            params.probe.name);
    end

    focus = [0, 0, params.elev_focus];

    Tx = xdc_linear_array(params.n_elements, params.element_width, ...
        params.element_height, params.kerf, ...
        params.n_sub_x, params.n_sub_y, focus);

    Rx = xdc_linear_array(params.n_elements, params.element_width, ...
        params.element_height, params.kerf, ...
        params.n_sub_x, params.n_sub_y, focus);

    t_ex = (0 : 1/params.fs : params.n_cycles/params.f0);
    excitation = sin(2 * pi * params.f0 * t_ex);
    xdc_excitation(Tx, excitation);

    t_ir = (0 : 1/params.fs : 2/params.f0);
    impulse = sin(2 * pi * params.f0 * t_ir) .* hanning(length(t_ir))';
    xdc_impulse(Tx, impulse);
    xdc_impulse(Rx, impulse);
end

function lri = das_beamform_planewave(rf_data, t_start, angle, params)
    [n_samples, n_elem] = size(rf_data);

    [X, Z] = meshgrid(params.x_grid, params.z_grid);
    elem_x = ((0:n_elem-1) - (n_elem-1)/2) * params.pitch;
    t_tx = (Z * cos(angle) + X * sin(angle)) / params.c;

    rf_padded = [rf_data; zeros(2, n_elem)];
    lri = zeros(size(X));

    for e = 1:n_elem
        t_rx = sqrt((X - elem_x(e)).^2 + Z.^2) / params.c;
        idx = (t_tx + t_rx - t_start) * params.fs + 1;

        idx_floor = floor(idx);
        frac = idx - idx_floor;

        valid = (idx_floor >= 1) & (idx_floor < n_samples);
        idx_safe = max(1, min(idx_floor, n_samples));

        val = (1 - frac) .* rf_padded(idx_safe + (e-1)*size(rf_padded,1)) ...
            + frac .* rf_padded(idx_safe + 1 + (e-1)*size(rf_padded,1));
        val(~valid) = 0;

        lri = lri + val;
    end
end

function bmode = envelope_to_bmode(envelope, dynamic_range)
    env_max = max(envelope(:));
    if env_max > 0
        bmode = 20 * log10(envelope / env_max);
    else
        bmode = zeros(size(envelope));
    end

    bmode = max(bmode, -dynamic_range);
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
    metadata.not_tracking_inputs = {'lri_frames', 'hri_frames', 'env_frames', 'bmode_frames', 'bmode_png', 'videos'};
end

function value = get_optional_param(params, name, default_value)
    if isfield(params, name)
        value = params.(name);
    else
        value = default_value;
    end
end

function save_lri_pngs(params, lri_bmode_frames, out_dir)
    lri_dir = fullfile(out_dir, 'lri_angle_png');
    if ~exist(lri_dir, 'dir'), mkdir(lri_dir); end

    for frame = 1:params.n_frames
        for a = 1:params.n_angles
            angle_deg = round(rad2deg(params.angles(a)));
            png_img = bmode_to_uint8(lri_bmode_frames(:,:,a,frame), params.dynamic_range);
            imwrite(png_img, fullfile(lri_dir, ...
                sprintf('frame_%03d_angle_%+03ddeg.png', frame, angle_deg)));
        end
    end

    save_lri_montage(params, lri_bmode_frames, out_dir);
end

function save_lri_montage(params, lri_bmode_frames, out_dir)
    figure('Name', 'Field II LRI angles', 'Color', 'w', ...
        'Position', [100, 100, 1200, 360]);
    for a = 1:params.n_angles
        angle_deg = round(rad2deg(params.angles(a)));
        subplot(1, params.n_angles, a);
        imagesc(params.x_grid * 1e3, params.z_grid * 1e3, lri_bmode_frames(:,:,a,1));
        colormap gray;
        axis image;
        xlabel('Lateral [mm]');
        ylabel('Depth [mm]');
        title(sprintf('%+d deg', angle_deg));
        clim([-params.dynamic_range 0]);
    end
    sgtitle('Field II uncompounded low-resolution images, frame 1');
    saveas(gcf, fullfile(out_dir, 'lri_angles_frame1.png'));
    close(gcf);
end

function img = bmode_to_uint8(bmode, dynamic_range)
    scaled = 255 * (bmode + dynamic_range) / dynamic_range;
    scaled(bmode <= -dynamic_range) = 0;
    scaled(bmode >= 0) = 255;
    img = uint8(scaled);
end
