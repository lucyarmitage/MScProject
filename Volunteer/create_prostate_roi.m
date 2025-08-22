clear;

S = load('C:/Users/lucya/MSC_PROJECT/Volunteer/recon_results/prostate_8mm_4001_short_skip100.mat');
Qmaps = S.Qmaps;

T1_recon = flip(fliplr(double(Qmaps(:,:,1,1))));
img = T1_recon;

% Display for ROI drawing
figure; imagesc(img,[0 1500]); axis image off; colorbar;

h = drawcircle('FaceAlpha',0.15,'Color','c','LineWidth',1.5);
wait(h);

% Create binary mask
roi_mask = createMask(h); 

hold on;
B = bwperim(roi_mask);
[y,x] = find(B);
plot(x,y,'c','LineWidth',1.5);
hold off;

save('prostate_roi_mask.mat','roi_mask');
disp('Circle ROI mask saved to roi_mask.mat');


