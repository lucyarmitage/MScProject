folder = 'C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/recon_results';
load(fullfile(folder, 'epg.mat'));

x = 77; 
y = 53;  
pixel_idx = sub2ind([N, N], y, x);

T1 = idx(idx2(pixel_idx), 1);
T2 = idx(idx2(pixel_idx), 2);
fprintf('Voxel (%d, %d): T1 = %.1f ms, T2 = %.1f ms\n', x, y, T1, T2);


% Measured signal
compressed_signal = squeeze(svd_images(y, x, :));
measured_signal = D.u * compressed_signal;

% Dictionary signal
dict_signal = dict0(:, idx2(pixel_idx));

% Scaling 
scaling_factor = norm(measured_signal) / norm(dict_signal);
dict_signal_scaled = dict_signal * scaling_factor;


figure;
plot(abs(measured_signal), 'k--', 'LineWidth', 1); hold on;
plot(abs(dict_signal_scaled), 'r', 'LineWidth', 1);
legend('Measured signal', 'Matched dictionary signal (scaled)');
xlabel('Timepoint'); ylabel('|Signal|');
title(sprintf('Abs - Voxel (%d, %d): T1=%d, T2=%d', x, y, T1, T2));

figure;
plot(real(measured_signal), 'k--', 'LineWidth', 1); hold on;
plot(real(dict_signal_scaled), 'r', 'LineWidth', 1);
legend('Measured signal', 'Matched dictionary signal (scaled)');
xlabel('Timepoint'); ylabel('|Signal|');
title(sprintf('Real - Voxel (%d, %d): T1=%d, T2=%d', x, y, T1, T2));

figure;
plot(imag(measured_signal), 'k--', 'LineWidth', 1); hold on;
plot(imag(dict_signal_scaled), 'r', 'LineWidth', 1);
legend('Measured signal', 'Matched dictionary signal (scaled)');
xlabel('Timepoint'); ylabel('|Signal|');
title(sprintf('Imag - Voxel (%d, %d): T1=%d, T2=%d', x, y, T1, T2));