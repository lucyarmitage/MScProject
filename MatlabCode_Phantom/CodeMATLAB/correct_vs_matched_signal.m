clear;

folder = 'C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/recon_results';
dict_files = {
    '10mm_101_short.mat'
    '10mm_301_short.mat'
    '10mm_501_short.mat'
    '10mm_701_short.mat'
};
dict_labels = {
    '101'
    '301'
    '501'
    '701'
};

x_flipped = 62; 
y_flipped = 51;

load(fullfile(folder, dict_files{1}));
x = N - x_flipped + 1;  
y = N - y_flipped + 1;  
pixel_idx = sub2ind([N, N], y, x);
compressed_meas_sig = squeeze(svd_images(y, x, :));
meas_sig = D.u * compressed_meas_sig;

matched_idx = idx2(pixel_idx);
T1_match = idx(matched_idx, 1);
T2_match = idx(matched_idx, 2);

gs_data = load('T1T2_results.mat');
t1_gs = imresize(gs_data.T1map, [N N]);
t2_gs = imresize(gs_data.T2map, [N N]);
gs_T1 = t1_gs(y_flipped, x_flipped);
gs_T2 = t2_gs(y_flipped, x_flipped);
fprintf('GS T1 = %.1f ms, GS T2 = %.1f ms\n', gs_T1, gs_T2);
fprintf('Matched T1 = %.1f ms, T2 = %.1f ms\n', T1_match, T2_match);

dists = sqrt((idx(:,1) - gs_T1).^2 + (idx(:,2) - gs_T2).^2);
[~, gs_idx] = min(dists);
gs_T1_dict = idx(gs_idx, 1);
gs_T2_dict = idx(gs_idx, 2);
fprintf('Closest dictionary entry to GS: T1 = %.1f, T2 = %.1f\n', gs_T1_dict, gs_T2_dict);

colors = lines(numel(dict_files));
figure_abs = figure('Name', 'Abs');
figure_real = figure('Name', 'Real');
figure_imag = figure('Name', 'Imag');

legend_entries_abs = {};
legend_entries_real = {};
legend_entries_imag = {};

for i = 1:numel(dict_files)
    file = dict_files{i};
    data = load(fullfile(folder, file));
    D = data.D;
    idx = data.idx;
    idx2 = data.idx2;
    dict0 = data.dict0;

    dict_sig = conj(dict0(:, matched_idx));
    scaling_factor = norm(meas_sig) / norm(dict_sig);
    dict_sig_scaled = dict_sig * scaling_factor;
    phi = angle(sum(conj(meas_sig) .* dict_sig_scaled));
    dict_sig_aligned = dict_sig_scaled * exp(-1i * phi);

    gs_dict_sig = dict0(:, gs_idx);
    scaling_factor_gs = norm(meas_sig) / norm(gs_dict_sig);
    gs_dict_sig_scaled = gs_dict_sig * scaling_factor_gs;
    phi_gs = angle(sum(conj(meas_sig) .* gs_dict_sig_scaled));
    gs_dict_sig_aligned = gs_dict_sig_scaled * exp(-1i * phi_gs);

    figure(figure_abs);
    plot(abs(dict_sig_aligned), '--', 'Color', colors(i,:), 'LineWidth', 1.2); hold on;
    plot(abs(gs_dict_sig_aligned), '-', 'Color', colors(i,:), 'LineWidth', 1.2);
    legend_entries_abs{end+1} = sprintf('%s Matched (T1=%.0f, T2=%.0f)', dict_labels{i}, T1_match, T2_match);
    legend_entries_abs{end+1} = [dict_labels{i} ' GS dict'];

    figure(figure_real);
    plot(real(dict_sig_aligned), '--', 'Color', colors(i,:), 'LineWidth', 1.2); hold on;
    plot(real(gs_dict_sig_aligned), '-', 'Color', colors(i,:), 'LineWidth', 1.2);
    legend_entries_real{end+1} = sprintf('%s Matched (T1=%.0f, T2=%.0f)', dict_labels{i}, T1_match, T2_match);
    legend_entries_real{end+1} = [dict_labels{i} ' GS dict'];

    figure(figure_imag);
    plot(imag(dict_sig_aligned), '--', 'Color', colors(i,:), 'LineWidth', 1.2); hold on;
    plot(imag(gs_dict_sig_aligned), '-', 'Color', colors(i,:), 'LineWidth', 1.2);
    legend_entries_imag{end+1} = sprintf('%s Matched (T1=%.0f, T2=%.0f)', dict_labels{i}, T1_match, T2_match);
    legend_entries_imag{end+1} = [dict_labels{i} ' GS dict'];
end

figure(figure_abs);
plot(abs(meas_sig), 'k', 'LineWidth', 2);
legend_entries_abs{end+1} = 'Measured';
xlabel('Timepoint'); ylabel('|Signal|');
title(sprintf('Abs: Voxel (%d,%d) | GS T1=%.0f, T2=%.0f | Closest dict: T1=%.0f, T2=%.0f', ...
    x, y, gs_T1, gs_T2, gs_T1_dict, gs_T2_dict));
legend(legend_entries_abs, 'Location', 'best'); grid on;

figure(figure_real);
plot(real(meas_sig), 'k', 'LineWidth', 2);
legend_entries_real{end+1} = 'Measured';
xlabel('Timepoint'); ylabel('Real');
title(sprintf('Real: Voxel (%d,%d) | GS T1=%.0f, T2=%.0f | Closest dict: T1=%.0f, T2=%.0f', ...
    x, y, gs_T1, gs_T2, gs_T1_dict, gs_T2_dict));
legend(legend_entries_real, 'Location', 'best'); grid on;

figure(figure_imag);
plot(imag(meas_sig), 'k', 'LineWidth', 2);
legend_entries_imag{end+1} = 'Measured';
xlabel('Timepoint'); ylabel('Imag');
title(sprintf('Imag: Voxel (%d,%d) | GS T1=%.0f, T2=%.0f | Closest dict: T1=%.0f, T2=%.0f', ...
    x, y, gs_T1, gs_T2, gs_T1_dict, gs_T2_dict));
legend(legend_entries_imag, 'Location', 'best'); grid on;
