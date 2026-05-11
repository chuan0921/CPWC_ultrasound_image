%% main_fieldii_parallel_simulation.m
% Field II run for the literature linear array.
% Parallelism is over angles within a seed (each worker owns its own Field II state).

clear; close all; clc;

proj_dir = fileparts(mfilename('fullpath'));
fieldii_dir = '/Users/dorishsu/MBP16/NTU/PCLAB/Project/m_files';
addpath(fieldii_dir);

params = setup_parameters('literature_l12_3v');
params.random_seeds = 1;
params.n_frames = 20;
params.max_parallel_workers = 5;  % 平行 5 angles
run_id = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
params.run_id = run_id;
params.run_timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));

if ~strcmp(params.probe.geometry, 'linear_array')
    error('This test script is intentionally limited to the literature linear array.');
end

out_root = fullfile(proj_dir, 'output_fieldii_parallel', params.probe.name, ['run_' run_id]);
if ~exist(out_root, 'dir'), mkdir(out_root); end
params.output_root = out_root;

fprintf('=== Parallel Field II Simulation ===\n');
fprintf('Probe: %s, seeds=%d, frames/seed=%d, angles=%d\n', ...
    params.probe.name, numel(params.random_seeds), params.n_frames, params.n_angles);
fprintf('Run ID: %s\n', params.run_id);
fprintf('Output: %s\n', out_root);

pool = gcp('nocreate');
if isempty(pool) || pool.NumWorkers < params.max_parallel_workers
    if ~isempty(pool); delete(pool); end
    try
        pool = parpool('Processes', params.max_parallel_workers);
    catch
        pool = parpool(params.max_parallel_workers);
    end
end
fprintf('Parallel pool: %d workers (parallelize over angles)\n', pool.NumWorkers);

seed_results = cell(numel(params.random_seeds), 1);
for seed_idx = 1:numel(params.random_seeds)
    seed_params = params;
    seed_params.random_seed = seed_params.random_seeds(seed_idx);
    seed_out_dir = fullfile(out_root, sprintf('seed_%04d', seed_params.random_seed));
    seed_results{seed_idx} = simulate_fieldii_seed(seed_params, seed_out_dir, fieldii_dir);
end

save(fullfile(out_root, 'batch_summary.mat'), 'seed_results', 'params');
fprintf('=== Field II Parallel Batch Complete ===\n');
