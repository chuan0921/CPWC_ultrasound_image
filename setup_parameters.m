function params = setup_parameters(probe_name)
% setup_parameters - 定義 CPWC 超音波流場模擬所有參數
% 參考: Leow & Tang, UMB 2018, Table 2

    if nargin < 1
        probe_name = 'literature_l12_3v';
    end

    %% 物理常數
    params.c = 1540;                        % 聲速 [m/s]

    %% 探頭參數
    % 可選: 'literature_l12_3v', 'linear_l12_5mhz' 或 'zipper_array'
    params.probe_name = probe_name;
    params = apply_probe_profile(params);

    %% 成像參數 (Compounded Plane Wave)
    params.f_transmit = params.f0;          % 發射頻率 [Hz]
    params.n_cycles = 1;                    % 激發脈衝週期數
    params.PRF = 10e3;                      % 脈衝重複頻率 [Hz]
    params.n_angles = 5;                    % Compounding 角度數
    params.angles = deg2rad([-10 -5 0 5 10]); % 平面波角度 [rad]

    %% 流場參數 (Poiseuille steady flow)
    params.phantom_name = 'literature_steady_poiseuille';
    params.v0 = 0.50;                       % 中心線速度 [m/s] (50 cm/s)
    params.R = 2.5e-3;                      % 管腔半徑 [m] (5 mm 內徑)
    params.wall_thickness = 0.2e-3;         % 管壁厚度 [m] (5.4 mm 外徑)
    params.vessel_center_z = 20e-3;         % 管中心深度 [m]
    params.beam_to_flow_angle = pi/2;       % 超音波束與流向夾角 [rad]

    %% 影像網格 (窄 FOV，只供 speckle tracking + 看血管壁)
    params.x_min = -5e-3;                   % FOV 寬 10 mm
    params.x_max =  5e-3;
    params.z_min = 16.3e-3;                 % 含外壁 z=17.3, 22.7 mm 各留 1 mm margin
    params.z_max = 23.7e-3;
    params.dx = params.lambda / 2;          % Lateral 像素間距 [m]
    params.dz = params.lambda / 4;          % Axial 像素間距 [m]

    params.x_grid = params.x_min : params.dx : params.x_max;
    params.z_grid = params.z_min : params.dz : params.z_max;
    params.Nx = length(params.x_grid);
    params.Nz = length(params.z_grid);

    %% 3D phantom elevation grid
    params.dy = params.lambda / 2;          % Elevation voxel spacing [m]
    params.elevation_beam_sigma = ...
        params.lambda * params.elev_focus / params.element_height; % Probe-determined elevation beam sigma [m]
    params.elevation_extent_half = 5e-3;    % Fixed linear-array 3D simulation support [m]
    if isfield(params.probe, 'element_data')
        y_vertices = params.probe.element_data(:, [3 6 9 12]);
        params.y_min = min(y_vertices(:));
        params.y_max = max(y_vertices(:));
    else
        params.y_min = -params.elevation_extent_half;
        params.y_max =  params.elevation_extent_half;
    end
    params.y_grid = params.y_min : params.dy : params.y_max;
    params.Ny = length(params.y_grid);

    %% 散射子參數
    % 管長度延伸到 FOV 外避免邊界效應
    params.vessel_length = 12e-3;           % 管長 [m] (wrap-around 安全下限: v0*T_total = 5 mm)
    params.scatter_density = 10;            % 每個解析度單元的散射子數
    params.wall_amp_factor = 10;            % 管壁相對振幅倍數

    %% 雜訊
    params.SNR_dB = 20;                     % 訊噪比 [dB]
    params.noise_mode = 'snr';              % 'snr' uses clean image RMS; 'floor' uses image_noise_floor
    params.bmode_reference_mode = 'seed_max'; % Fixed dB reference per seed for temporal consistency

    %% 模擬設定
    params.random_seeds = 1:5;              % 批次產生不同 speckle seed
    params.random_seed = params.random_seeds(1);  % 單次模擬預設 seed
    params.image_noise_floor = 0.01;        % 影像域複數雜訊
    params.n_frames = 20;                   % 每個 seed 的模擬幀數
    params.video_fps = 20;                  % 輸出影片幀率 [frames/s]
    params.dynamic_range = 50;              % B-mode 動態範圍 [dB]

end

function params = apply_probe_profile(params)
    probe = probe_profile(params.probe_name);
    params.probe = probe;

    params.f0 = probe.f0;
    params.lambda = params.c / params.f0;
    params.n_elements = probe.n_elements;
    params.pitch = probe.pitch;
    params.kerf = probe.kerf;
    params.element_width = probe.element_width;
    params.element_height = probe.element_height;
    params.elev_focus = probe.elev_focus;
    params.n_sub_x = probe.n_sub_x;
    params.n_sub_y = probe.n_sub_y;
    params.fs = probe.fs;

    if isfield(probe, 'element_data')
        params.element_data = probe.element_data;
    end
    if isfield(probe, 'element_centers')
        params.element_centers = probe.element_centers;
    end
end

function probe = probe_profile(name)
    key = lower(strrep(name, '-', '_'));

    switch key
        case {'literature_l12_3v', 'l12_3v', 'paper'}
            probe.name = 'literature_l12_3v';
            probe.label = 'Literature L12-3v';
            probe.geometry = 'linear_array';
            probe.f0 = 8e6;
            probe.n_elements = 128;
            probe.pitch = 0.2e-3;
            probe.kerf = 0.02e-3;
            probe.element_height = 5e-3;
            probe.elev_focus = 20e-3;
            probe.n_sub_x = 1;
            probe.n_sub_y = 1;
            probe.fs = 100e6;

        case {'linear_l12_5mhz', 'linear_5mhz', 'l12_5mhz'}
            probe.name = 'linear_l12_5mhz';
            probe.label = 'Linear L12 geometry, 5 MHz';
            probe.geometry = 'linear_array';
            probe.f0 = 5e6;
            probe.n_elements = 128;
            probe.pitch = 0.2e-3;
            probe.kerf = 0.02e-3;
            probe.element_height = 5e-3;
            probe.elev_focus = 20e-3;
            probe.n_sub_x = 1;
            probe.n_sub_y = 1;
            probe.fs = 100e6;

        case {'zipper_array', 'zipper'}
            probe.name = 'zipper_array';
            probe.label = 'Zipper array';
            probe.geometry = 'zipper_array';
            probe.f0 = 5e6;
            probe.gap = 0.01e-3;
            probe.n_x = 128;
            probe.n_y = 2;
            probe.n_elements = probe.n_x * probe.n_y;
            probe.pitch = 0.31e-3;
            probe.kerf = probe.gap;
            probe.element_width = 0.3e-3;
            probe.element_height = 20e-3;
            probe.row_pitch = 20e-3 + probe.gap;
            probe.elev_focus = NaN;
            probe.n_sub_x = 1;
            probe.n_sub_y = 1;
            probe.fs = 100e6;

            probe.element_data = zipper_element_data(probe);
            probe.element_centers = element_centers_from_data(probe.element_data);

        otherwise
            error('Unknown probe profile: %s', name);
    end

    if ~isfield(probe, 'element_width')
        probe.element_width = probe.pitch - probe.kerf;
    end
end

function data = zipper_element_data(probe)
    n = probe.n_elements;
    data = zeros(n, 19);

    aperture_width = (probe.n_x - 1) * probe.pitch + probe.element_width;
    x0 = -aperture_width / 2;
    y_centers = [-0.5, 0.5] * probe.row_pitch;

    idx = 0;
    for ix = 1:probe.n_x
        x_left = x0 + (ix - 1) * probe.pitch;
        x_right = x_left + probe.element_width;

        for iy = 1:probe.n_y
            idx = idx + 1;
            y_center = y_centers(iy);
            y_bottom = y_center - probe.element_height / 2;
            y_top = y_center + probe.element_height / 2;

            data(idx,:) = [idx, ...
                x_left,  y_top,    0, ...
                x_left,  y_bottom, 0, ...
                x_right, y_bottom, 0, ...
                x_right, y_top,    0, ...
                1, probe.element_width, probe.element_height, 0, 0, 0];
        end
    end
end

function centers = element_centers_from_data(data)
    centers = zeros(size(data, 1), 3);
    for i = 1:size(data, 1)
        x_coords = [data(i,2), data(i,5), data(i,8), data(i,11)];
        y_coords = [data(i,3), data(i,6), data(i,9), data(i,12)];
        z_coords = [data(i,4), data(i,7), data(i,10), data(i,13)];
        centers(i,:) = [mean(x_coords), mean(y_coords), mean(z_coords)];
    end
end
