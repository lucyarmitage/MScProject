ref = load('Prostate_Ref_T1_T2.mat');
recon = load('C:/Users/lucya/MSC_PROJECT/Volunteer/recon_results/prostate_20mm_101_short.mat');

T1_ref = ref.T1map;
T2_ref = ref.T2map;

Qmaps = recon.Qmaps; 
T1_recon = flip(fliplr(Qmaps(:,:,1,1)));
T2_recon = flip(fliplr(Qmaps(:,:,1,2)));


figure;
subplot(2,3,1); imagesc(T1_ref); axis image off; colorbar; title('T1 reference');
subplot(2,3,2); imagesc(T1_recon); axis image off; colorbar; title('T1 reconstructed');


load('prostate_square.mat');

T1_recon_roi = T1_recon(r1:r2, c1:c2);
T2_recon_roi = T2_recon(r1:r2, c1:c2);

T1_recon_mean = mean(T1_recon_roi(:));
T2_recon_mean = mean(T2_recon_roi(:));

fprintf('T1 = %.2f ms\n', T1_recon_mean);
fprintf('T2 = %.2f ms\n', T2_recon_mean);