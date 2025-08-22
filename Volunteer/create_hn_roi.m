clear;

S = load('C:/Users/lucya/MSC_PROJECT/Volunteer/recon_results/HN_20mm_101_short.mat');
Qmaps = S.Qmaps;

T1_recon = flip(fliplr(double(Qmaps(:,:,1,1))));
img = T1_recon;

figure; imagesc(img,[0 1500]); axis image off;

h1 = drawcircle('FaceAlpha',0,'Color','w','LineWidth',1.5);
wait(h1);
mask1 = createMask(h1);

h2 = drawcircle('FaceAlpha',0,'Color','w','LineWidth',1.5); 
wait(h2);
mask2 = createMask(h2);

% Combine 2 masks 
roi_mask = mask1 | mask2;

hold on;
B = bwperim(roi_mask_HN);
[y,x] = find(B);
plot(x,y,'w','LineWidth',1.5);
hold off;

save('roi_mask_HN.mat','roi_mask');
