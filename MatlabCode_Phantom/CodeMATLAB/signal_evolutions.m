folder = 'C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/recon_results';
load(fullfile(folder, '20mm_101.mat'));

x = 77; 
y = 53;  
pixel_idx = sub2ind([N, N], y, x);

T1 = idx(idx2(pixel_idx), 1);
T2 = idx(idx2(pixel_idx), 2);
fprintf('Voxel (%d, %d): T1 = %.1f ms, T2 = %.1f ms\n', x, y, T1, T2);

% Measured signal
compressed_meas_sig = squeeze(svd_images(y, x, :));
meas_sig = D.u * compressed_meas_sig;

% Dictionary signal
compressed_dict_sig = D.magnetization(:, idx2(pixel_idx));
dict_sig = dict0(:, idx2(pixel_idx));

% Scaling 
scaling_factor_compr = norm(compressed_meas_sig) / norm(compressed_dict_sig);
compressed_dict_sig_scaled = compressed_dict_sig * scaling_factor_compr;

scaling_factor = norm(meas_sig) / norm(dict_sig);
dict_sig_scaled = dict_sig * scaling_factor;



% Compressed signal

figure;
plot(abs(compressed_meas_sig), 'k--'); hold on;
plot(abs(compressed_dict_sig_scaled), 'r');
legend('Measured signal', "Matched dictionary signal (scaled)")
xlabel('Timepoint'); ylabel('|Signal|');                                       % Fix axis labels !!
title(sprintf('Abs - Voxel (%d, %d): T1=%d, T2=%d', x, y, T1, T2));

figure;
plot(real(compressed_meas_sig), 'k--'); hold on;
plot(real(compressed_dict_sig_scaled), 'r');
legend('Measured signal', "Matched dictionary signal (scaled)")
xlabel('Timepoint'); ylabel('|Signal|');                                       % Fix axis labels !!
title(sprintf('Real - Voxel (%d, %d): T1=%d, T2=%d', x, y, T1, T2));

figure;
plot(imag(compressed_meas_sig), 'k--'); hold on;
plot(imag(compressed_dict_sig_scaled), 'r');
legend('Measured signal', "Matched dictionary signal (scaled)")
xlabel('Timepoint'); ylabel('|Signal|');                                       % Fix axis labels !!
title(sprintf('Imag - Voxel (%d, %d): T1=%d, T2=%d', x, y, T1, T2));

% Uncompressed signal 

figure;
plot(abs(meas_sig), 'k--'); hold on;
plot(abs(dict_sig_scaled), 'r');
legend('Measured signal', 'Matched dictionary signal (scaled)');
xlabel('Timepoint'); ylabel('|Signal|');
title(sprintf('Abs - Voxel (%d, %d): T1=%d, T2=%d', x, y, T1, T2));

figure;
plot(real(meas_sig), 'k--'); hold on;
plot(real(dict_sig_scaled), 'r');
legend('Measured signal', 'Matched dictionary signal (scaled)');
xlabel('Timepoint'); ylabel('|Signal|');
title(sprintf('Real - Voxel (%d, %d): T1=%d, T2=%d', x, y, T1, T2));

figure;
plot(imag(meas_sig), 'k--'); hold on;
plot(imag(dict_sig_scaled), 'r');
legend('Measured signal', 'Matched dictionary signal (scaled)');
xlabel('Timepoint'); ylabel('|Signal|');
title(sprintf('Imag - Voxel (%d, %d): T1=%d, T2=%d', x, y, T1, T2));

