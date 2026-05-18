%% plot_probe_geometry.m
% Plot the current receive/transmit linear-array probe geometry.

clear; close all; clc;

params = setup_parameters('zipper_array');

proj_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(proj_dir, 'output', params.probe.name);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

if isfield(params, 'element_centers')
    elem_x = params.element_centers(:,1).';
else
    elem_x = ((0:params.n_elements-1) - (params.n_elements-1)/2) * params.pitch;
end
elem_width = params.element_width;
elem_height = params.element_height;
aperture_width = max(elem_x) - min(elem_x) + elem_width;

fprintf('=== Probe Geometry ===\n');
fprintf('Probe profile: %s (%s)\n', params.probe.name, params.probe.label);
fprintf('Type: linear array\n');
fprintf('Elements: %d\n', params.n_elements);
fprintf('Pitch: %.3f mm\n', params.pitch * 1e3);
fprintf('Element width: %.3f mm\n', elem_width * 1e3);
fprintf('Kerf: %.3f mm\n', params.kerf * 1e3);
fprintf('Element height: %.3f mm\n', elem_height * 1e3);
fprintf('Physical aperture width: %.3f mm\n', aperture_width * 1e3);
if isfinite(params.elev_focus)
    fprintf('Elevation lens focus: %.3f mm\n', params.elev_focus * 1e3);
else
    fprintf('Elevation lens focus: not specified for this profile\n');
end

figure('Name', 'Zipper Array XY Geometry', 'Color', 'w', ...
    'Position', [100, 100, 900, 900]);

hold on;
if isfield(params, 'element_data')
    for e = 1:params.n_elements
        x_coords = [params.element_data(e,2), params.element_data(e,5), ...
            params.element_data(e,8), params.element_data(e,11)] * 1e3;
        y_coords = [params.element_data(e,3), params.element_data(e,6), ...
            params.element_data(e,9), params.element_data(e,12)] * 1e3;
        patch(x_coords, y_coords, [0.20 0.44 0.78], ...
            'EdgeColor', [0.12 0.24 0.42]);
        if params.n_elements <= 16 || e <= 4 || e > params.n_elements - 4
            text(mean(x_coords), mean(y_coords), sprintf('%d', e), ...
                'HorizontalAlignment', 'center', 'Color', 'w', ...
                'FontWeight', 'bold', 'FontSize', 7);
        end
    end
else
    for e = 1:params.n_elements
        x_left = elem_x(e) - elem_width / 2;
        rectangle('Position', [x_left * 1e3, -0.5, elem_width * 1e3, 1], ...
            'FaceColor', [0.20 0.44 0.78], 'EdgeColor', 'none');
    end
end
hold off;
axis image;
xlabel('Lateral x [mm]');
ylabel('Elevation y [mm]');
title(sprintf('%s x-y aperture: %d elements, f0 %.1f MHz', ...
    params.probe.label, params.n_elements, params.f0 / 1e6));
grid on;

filename = fullfile(out_dir, 'zipper_array_xy_elements.png');
saveas(gcf, filename);
fprintf('Saved: %s\n', filename);
