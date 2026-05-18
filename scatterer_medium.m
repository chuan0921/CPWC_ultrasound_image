function varargout = scatterer_medium(action, varargin)
% scatterer_medium - Probe-independent scatterer phantom utilities.
% Usage:
%   [flow_pos, flow_amp, static_pos, static_amp] = scatterer_medium('generate', params)
%   flow_pos = scatterer_medium('update_flow', flow_pos, params, dt)

    switch lower(action)
        case 'generate'
            [varargout{1:nargout}] = generate_scatterers(varargin{:});
        case 'update_flow'
            varargout{1} = update_flow_positions(varargin{:});
        case 'ground_truth'
            varargout{1} = generate_ground_truth(varargin{:});
        otherwise
            error('Unknown scatterer_medium action: %s', action);
    end
end

function [flow_pos, flow_amp, static_pos, static_amp] = generate_scatterers(params)
    R = params.R;
    z_c = params.vessel_center_z;

    F_elev = params.elev_focus / params.element_height;
    y_half = params.lambda * F_elev * 3;

    delta_x = params.lambda * 2;
    delta_z = params.lambda;
    delta_y = params.lambda * F_elev;
    V_cell = delta_x * delta_z * delta_y;

    y_clip = min(y_half, R);
    V_flow_full = pi * R^2 * params.vessel_length;
    N_flow_full = round(params.scatter_density * V_flow_full / V_cell);

    oversample = ceil(R / y_clip) + 1;
    N_gen = N_flow_full * oversample;
    r = R * sqrt(rand(N_gen, 1));
    theta = 2 * pi * rand(N_gen, 1);
    fy = r .* cos(theta);
    fz = r .* sin(theta);

    keep = abs(fy) <= y_half;
    fy = fy(keep);
    fz = fz(keep);
    N_flow = min(length(fy), N_flow_full);
    fy = fy(1:N_flow);
    fz = fz(1:N_flow);

    flow_x = (rand(N_flow, 1) - 0.5) * params.vessel_length;
    flow_pos = [flow_x, fy, z_c + fz];
    flow_amp = randn(N_flow, 1);

    wt = params.wall_thickness;
    R_outer = R + wt;
    V_wall_full = pi * (R_outer^2 - R^2) * params.vessel_length;
    N_wall_full = round(params.scatter_density * V_wall_full / V_cell);

    N_gen_w = N_wall_full * oversample;
    r_wall = sqrt(R^2 + (R_outer^2 - R^2) * rand(N_gen_w, 1));
    theta_wall = 2 * pi * rand(N_gen_w, 1);
    wy = r_wall .* cos(theta_wall);
    wz = r_wall .* sin(theta_wall);

    keep_w = abs(wy) <= y_half;
    wy = wy(keep_w);
    wz = wz(keep_w);
    N_wall = min(length(wy), N_wall_full);
    wy = wy(1:N_wall);
    wz = wz(1:N_wall);

    wall_x = (rand(N_wall, 1) - 0.5) * params.vessel_length;
    wall_pos = [wall_x, wy, z_c + wz];
    wall_amp = params.wall_amp_factor * randn(N_wall, 1);

    x_range = params.x_max - params.x_min;
    z_range = params.z_max - params.z_min;
    y_range = 2 * y_half;
    V_fov = x_range * y_range * z_range;
    V_vessel_in_fov = pi * R_outer^2 * x_range;
    V_tissue = V_fov - V_vessel_in_fov;
    N_tissue = round(params.scatter_density * V_tissue / V_cell);

    N_gen_t = round(N_tissue * 1.3);
    tissue_x = params.x_min + x_range * rand(N_gen_t, 1);
    tissue_y = (rand(N_gen_t, 1) - 0.5) * y_range;
    tissue_z = params.z_min + z_range * rand(N_gen_t, 1);

    dist_to_axis = sqrt(tissue_y.^2 + (tissue_z - z_c).^2);
    outside = dist_to_axis > R_outer;
    tissue_x = tissue_x(outside);
    tissue_y = tissue_y(outside);
    tissue_z = tissue_z(outside);

    if length(tissue_x) > N_tissue
        tissue_x = tissue_x(1:N_tissue);
        tissue_y = tissue_y(1:N_tissue);
        tissue_z = tissue_z(1:N_tissue);
    end

    tissue_pos = [tissue_x, tissue_y, tissue_z];
    tissue_amp = randn(length(tissue_x), 1);

    static_pos = [wall_pos; tissue_pos];
    static_amp = [wall_amp; tissue_amp];

    fprintf('Scatterers: flow=%d, wall=%d, tissue=%d, total=%d (y_half=%.1fmm)\n', ...
        N_flow, N_wall, length(tissue_x), N_flow + size(static_pos, 1), y_half*1e3);
end

function flow_pos = update_flow_positions(flow_pos, params, dt)
    R = params.R;
    z_c = params.vessel_center_z;
    half_len = params.vessel_length / 2;

    r = sqrt(flow_pos(:,2).^2 + (flow_pos(:,3) - z_c).^2);
    vx = params.v0 * max(0, 1 - (r / R).^2);

    flow_pos(:,1) = flow_pos(:,1) + vx * dt;

    out_right = flow_pos(:,1) > half_len;
    flow_pos(out_right, 1) = flow_pos(out_right, 1) - params.vessel_length;

    out_left = flow_pos(:,1) < -half_len;
    flow_pos(out_left, 1) = flow_pos(out_left, 1) + params.vessel_length;
end

function gt = generate_ground_truth(params)
    row_centers = row_y_centers(params);
    if numel(row_centers) >= 2
        row_separation = abs(row_centers(2) - row_centers(1));
    else
        row_separation = NaN;
    end

    gt.vessel.center_x = params.vessel_center_x;
    gt.vessel.center_z = params.vessel_center_z;
    gt.vessel.radius = params.R;
    gt.vessel.diameter = params.R * 2;
    gt.vessel.wall_upper = params.vessel_center_z - params.R;
    gt.vessel.wall_lower = params.vessel_center_z + params.R;
    gt.vessel.wall_thickness = params.wall_thickness;

    gt.flow.v0 = params.v0;
    gt.flow.axis = params.flow_axis;
    gt.flow.beam_to_flow_angle = params.beam_to_flow_angle;
    gt.flow.viscosity = 0.004;
    gt.flow.WSR_analytical = 2 * params.v0 / params.R;
    gt.flow.Re = params.v0 * params.R * 2 * 1060 / gt.flow.viscosity;
    gt.flow.Q_truth = pi * params.R^2 * params.v0 / 2;
    gt.flow.row_y_centers = row_centers;
    gt.flow.row_separation = row_separation;

    [X, Z] = meshgrid(params.x_grid, params.z_grid);
    r_xz = sqrt((X - params.vessel_center_x).^2 + ...
        (Z - params.vessel_center_z).^2);
    inside_xz = r_xz <= params.R;

    vy_map_xz = zeros(size(X));
    vy_map_xz(inside_xz) = params.v0 * ...
        (1 - (r_xz(inside_xz) / params.R).^2);

    gt.flow.vy_map_xz = vy_map_xz;
    gt.flow.vessel_mask_xz = inside_xz;

    [~, ix_center] = min(abs(params.x_grid));
    gt.flow.radial_profile_z = params.z_grid';
    gt.flow.radial_profile_v = vy_map_xz(:, ix_center);

    dt_frame = params.n_angles / params.PRF;
    gt.flow.dt_frame = dt_frame;
    gt.flow.dy_map_xz = vy_map_xz * dt_frame;
    gt.flow.dy_max = params.v0 * dt_frame;
    gt.flow.expected_lag_s = row_separation / params.v0;
    gt.flow.expected_lag_frames = gt.flow.expected_lag_s / dt_frame;
    gt.flow.expected_lag_pulses = gt.flow.expected_lag_s * params.PRF;

    gt.grid.x = params.x_grid;
    gt.grid.y = params.y_grid;
    gt.grid.z = params.z_grid;
    gt.grid.dx = params.dx;
    gt.grid.dy = params.dy;
    gt.grid.dz = params.dz;
    gt.grid.Nx = params.Nx;
    gt.grid.Ny = params.Ny;
    gt.grid.Nz = params.Nz;

    fprintf(['Ground truth: vy0=%.2f m/s, R=%.1f mm, Q=%.3f mL/s, ', ...
        'row_sep=%.2f mm, lag=%.1f frames\n'], ...
        params.v0, params.R*1e3, gt.flow.Q_truth*1e6, ...
        row_separation*1e3, gt.flow.expected_lag_frames);
end

function centers = row_y_centers(params)
    if ~isfield(params.probe, 'n_y') || params.probe.n_y < 1 || ...
            ~isfield(params.probe, 'element_data')
        centers = NaN;
        return;
    end

    centers = zeros(1, params.probe.n_y);
    for row = 1:params.probe.n_y
        element_rows = row:params.probe.n_y:params.probe.n_elements;
        y_vertices = params.probe.element_data(element_rows, [3 6 9 12]);
        centers(row) = mean([min(y_vertices(:)), max(y_vertices(:))]);
    end
end
