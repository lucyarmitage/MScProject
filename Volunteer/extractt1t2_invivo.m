clear;

S = load('C:/Users/lucya/MSC_PROJECT/Volunteer/recon_results/prostate_8mm_2001_long.mat'); 
Qmaps = S.Qmaps;

T1_recon = flip(fliplr(double(Qmaps(:,:,1,1))));
T2_recon = flip(fliplr(double(Qmaps(:,:,1,2))));

M = load('prostate_roi_mask.mat'); % prostate
% M = load('roi_mask_HN.mat');   % HN
roi_mask = M.roi_mask;

% Extract values inside ROI
T1_vals = T1_recon(roi_mask);
T2_vals = T2_recon(roi_mask);

mean_T1 = mean(T1_vals, 'omitnan');
mean_T2 = mean(T2_vals, 'omitnan');
std_T1  = std(T1_vals,  'omitnan');
std_T2  = std(T2_vals,  'omitnan');

fprintf('Mean T1 in ROI: %.2f ms (± %.2f)\n', mean_T1, std_T1);
fprintf('Mean T2 in ROI: %.2f ms (± %.2f)\n', mean_T2, std_T2);

figure;
imagesc(T1_recon, [0 1500]); axis image off; colorbar;
title('T1 with circle ROI overlay');
hold on;
B = bwperim(roi_mask);
[y,x] = find(B);
plot(x,y,'c','LineWidth',1.5);
hold off;
