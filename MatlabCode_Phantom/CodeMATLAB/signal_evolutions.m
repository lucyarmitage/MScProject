% Load precomputed reconstruction and dictionary 
folder = 'C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/';
load(fullfile(folder, 'recon_results_epg.mat'));

% Choose voxel
x = 77; 
y = 53;  
pixel_idx = sub2ind([N, N], y, x);

T1 = idx(idx2(pixel_idx), 1);
T2 = idx(idx2(pixel_idx), 2);
fprintf('Voxel (%d, %d): T1 = %.1f ms, T2 = %.1f ms\n', x, y, T1, T2);

% Extract measured signal
compressed_signal = squeeze(svd_images(y, x, :));
measured_signal = D.u * compressed_signal;

% Matched dictionary signal
dict_signal = dict0(:, idx2(pixel_idx));
matched_T1 = idx(idx2(pixel_idx), 1);
matched_T2 = idx(idx2(pixel_idx), 2);

% Scale dictionary signal to match measured signal amplitude
scaling_factor = (measured_signal' * dict_signal) / (dict_signal' * dict_signal);
dict_signal_scaled = dict_signal * scaling_factor;


figure;
plot(abs(measured_signal), 'b', 'LineWidth', 1.5); hold on;
plot(abs(dict_signal_scaled), 'r--', 'LineWidth', 1.5);
legend('Measured signal', 'Matched dictionary signal (scaled)');
xlabel('Timepoint'); ylabel('|Signal|');
title(sprintf('Voxel (%d, %d)', x, y));

figure;
plot(real(measured_signal), 'b', 'LineWidth', 1.5); hold on;
plot(real(dict_signal_scaled), 'r--', 'LineWidth', 1.5);
legend('Measured signal', 'Matched dictionary signal (scaled)');
xlabel('Timepoint'); ylabel('|Signal|');
title(sprintf('Voxel (%d, %d)', x, y));

figure;
plot(imag(measured_signal), 'b', 'LineWidth', 1.5); hold on;
plot(imag(dict_signal_scaled), 'r--', 'LineWidth', 1.5);
legend('Measured signal', 'Matched dictionary signal (scaled)');
xlabel('Timepoint'); ylabel('|Signal|');
title(sprintf('Voxel (%d, %d)', x, y));
