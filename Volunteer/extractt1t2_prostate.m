ref = load('Prostate_Ref_T1_T2.mat');
recon = load('C:/Users/lucya/MSC_PROJECT/Volunteer/recon_results/prostate_20mm_101_short_zoom.mat');

T1_ref = ref.T1map;
T2_ref = ref.T2map;

Qmaps = recon.Qmaps; 
T1_recon = flip(fliplr(Qmaps(:,:,1,1)));
T2_recon = flip(fliplr(Qmaps(:,:,1,2)));
T1_recon_resized = imresize(T1_recon, size(T1_ref));
T2_recon_resized = imresize(T2_recon, size(T2_ref));

T1_diff = T1_recon_resized - T1_ref;
T2_diff = T2_recon_resized - T2_ref;

T1_abs_error = mean(abs(T1_diff(:)));
T1_rmse = sqrt(mean((T1_diff(:)).^2));
T2_abs_error = mean(abs(T2_diff(:)));
T2_rmse = sqrt(mean((T2_diff(:)).^2));

fprintf('T1 mean absolute error: %.2f ms\n', T1_abs_error);
fprintf('T1 RMSE: %.2f ms\n', T1_rmse);
fprintf('T2 mean absolute error: %.2f ms\n', T2_abs_error);
fprintf('T2 RMSE: %.2f ms\n', T2_rmse);

figure;
subplot(2,3,1); imagesc(T1_ref); axis image off; colorbar; title('T1 reference');
subplot(2,3,2); imagesc(T1_recon_resized); axis image off; colorbar; title('T1 reconstructed');
subplot(2,3,3); imagesc(T1_diff); axis image off; colorbar; title('T1 error map');

subplot(2,3,4); imagesc(T2_ref); axis image off; colorbar; title('T2 reference');
subplot(2,3,5); imagesc(T2_recon_resized); axis image off; colorbar; title('T2 reconstructed');
subplot(2,3,6); imagesc(T2_diff); axis image off; colorbar; title('T2 error map');


load('prostate_square.mat');

T1_recon_roi = T1_recon_resized(r1:r2, c1:c2);
T2_recon_roi = T2_recon_resized(r1:r2, c1:c2);

T1_recon_mean = mean(T1_recon_roi(:));
T2_recon_mean = mean(T2_recon_roi(:));

fprintf('T1 = %.2f ms\n', T1_recon_mean);
fprintf('T2 = %.2f ms\n', T2_recon_mean);