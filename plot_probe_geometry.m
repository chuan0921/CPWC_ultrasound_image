%% plot_probe_geometry.m
% Plot the current receive/transmit linear-array probe geometry.

clear; close all; clc;

params = setup_parameters();

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

figure('Name', 'Probe Geometry', 'Color', 'w', 'Position', [100, 100, 1250, 560]);

subplot(2, 1, 1);
hold on;
if isfield(params, 'element_data')
    for e = 1:params.n_elements
        x_coords = [params.element_data(e,2), params.element_data(e,5), ...
            params.element_data(e,8), params.element_data(e,11)] * 1e3;
        y_coords = [params.element_data(e,3), params.element_data(e,6), ...
            params.element_data(e,9), params.element_data(e,12)] * 1e3;
        patch(x_coords, y_coords, [0.20 0.44 0.78], ...
            'EdgeColor', [0.12 0.24 0.42]);
        text(mean(x_coords), mean(y_coords), sprintf('%d', e), ...
            'HorizontalAlignment', 'center', 'Color', 'w', 'FontWeight', 'bold');
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
title(sprintf('%s: %d elements, f0 %.1f MHz', ...
    params.probe.label, params.n_elements, params.f0 / 1e6));
grid on;

subplot(2, 2, 3);
hold on;
if isfield(params, 'element_data')
    for e = 1:params.n_elements
        x_coords = [params.element_data(e,2), params.element_data(e,5), ...
            params.element_data(e,8), params.element_data(e,11)] * 1e3;
        y_coords = [params.element_data(e,3), params.element_data(e,6), ...
            params.element_data(e,9), params.element_data(e,12)] * 1e3;
        patch(x_coords, y_coords, [0.20 0.44 0.78], ...
            'EdgeColor', [0.12 0.24 0.42]);
    end
else
    for e = 1:params.n_elements
        x_left = elem_x(e) - elem_width / 2;
        rectangle('Position', [x_left * 1e3, -elem_height * 0.5e3, ...
            elem_width * 1e3, elem_height * 1e3], ...
            'FaceColor', [0.20 0.44 0.78], 'EdgeColor', [0.12 0.24 0.42]);
    end
end
hold off;
axis image;
xlabel('Lateral x [mm]');
ylabel('Elevation y [mm]');
title('Aperture face: x-y view');
grid on;

subplot(2, 2, 4);
hold on;
plot(elem_x * 1e3, zeros(size(elem_x)), '.', 'MarkerSize', 8, ...
    'DisplayName', 'elements');

z_lines = linspace(params.z_min, params.z_max, 4);
if z_lines(1) == 0
    z_lines(1) = params.z_max / 8;
end
x_span = [-aperture_width aperture_width] * 0.55;
colors = lines(length(params.angles));
for a = 1:length(params.angles)
    angle = params.angles(a);
    for zi = 1:length(z_lines)
        z0 = z_lines(zi);
        x_center = z0 * tan(angle);
        plot((x_center + x_span) * 1e3, [z0 z0] * 1e3, ...
            'Color', colors(a,:), 'LineWidth', 0.9);
    end
end
hold off;
axis image;
xlabel('Lateral x [mm]');
ylabel('Depth z [mm]');
title('Plane-wave transmit fronts in x-z (no lateral focus)');
grid on;

annotation('textbox', [0.54 0.03 0.42 0.08], ...
    'String', plane_wave_note(params), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left');

saveas(gcf, fullfile(out_dir, 'probe_geometry.png'));
fprintf('Saved: %s\n', fullfile(out_dir, 'probe_geometry.png'));

function note = plane_wave_note(params)
    if isfinite(params.elev_focus)
        note = sprintf(['Plane-wave angles are steered by element delays; ', ...
            'there is no x-z transmit focus. The %.1f mm focus is the fixed ', ...
            'elevation lens focus used by the linear array model.'], params.elev_focus * 1e3);
    else
        note = ['Plane-wave angles are steered by element delays; ', ...
            'there is no x-z transmit focus. This probe profile uses explicit ', ...
            'element polygons instead of a linear-array elevation focus.'];
    end
end
