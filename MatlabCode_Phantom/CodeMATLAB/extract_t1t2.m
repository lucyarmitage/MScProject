% --- Load precomputed reconstruction and dictionary ---
folder = 'C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/';
addpath(genpath('C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/CodeMATLAB/utils'));
load(fullfile(folder, 'recon_results101.mat'));
load('vial_mask_geometry.mat');  % loads 'centres' and 'radii'

% Qmaps = reshape(D.lookup_table(idx2,:),[[N N 1], size(D.lookup_table,2)]);
% Qmaps = cat(numel(size(Qmaps)),Qmaps,reshape(c ./ D.normalization(idx2).',[N N]));
% % PD    = c ./ D.normalization(idx2).';
% % PD    = reshape(PD, [N N]);
% 
% % Do some background masking
% mask = sum(abs(reshape(svd_images,[[N N] R])),3);
% mask(mask < .4* max(mask(:))) = 0;
% mask (mask ~= 0) = 1;
% Qmaps = bsxfun(@times,Qmaps,mask);
% 
% Qmaps = flip (Qmaps);
% % figure; imshow3(fliplr(Qmaps(:,:,1,1)), [0 max(max(Qmaps(:,:,1,1)))])
% % figure; imshow3(fliplr(Qmaps(:,:,1,2)), [0 max(max(Qmaps(:,:,1,2)))])
% % max(max(Qmaps(:,:,1,1)))
% % max(max(Qmaps(:,:,1,2)))
% 
% figure; imshow3(fliplr(Qmaps(:,:,1,1)), [0 3000]);
% figure; imshow3(fliplr(Qmaps(:,:,1,2)), [0 1500]); 
% % figure; imshow3(fliplr(Qmaps(:,:,1,1)), [0 3000]); colormap hot;
% % figure; imshow3(fliplr(Qmaps(:,:,1,2)), [0 1500]); colormap turbo;

magImg = fliplr(Qmaps(:,:,1,1))
  % Flip vertically to match Qmaps
[H, W] = size(magImg);
[xx, yy] = meshgrid(1:W, 1:H);

nVials = size(centres, 1);
Masks = false(H, W, nVials);

for k = 1:nVials
    cx = centres(k,1);
    cy = centres(k,2);
    r  = radii(k);
    Masks(:,:,k) = ((xx - cx).^2 + (yy - cy).^2) <= r^2;
end



T1 = fliplr(Qmaps(:,:,1,1));
T2 = fliplr(Qmaps(:,:,1,2));

T1_mean = zeros(nVials,1);
T2_mean = zeros(nVials,1);

for k = 1:nVials
    voxels = find(Masks(:,:,k));
    T1_mean(k) = mean(T1(voxels));
    T2_mean(k) = mean(T2(voxels));
end

results = table((1:nVials)', T1_mean, T2_mean, ...
    'VariableNames', {'Vial','T1_mean','T2_mean'});

disp(results)


figure; imshow(magImg, []); title('Final vial positions');
hold on;
viscircles(centres, radii, 'Color','y');
for k = 1:nVials
    text(centres(k,1), centres(k,2), sprintf('%d',k), ...
         'Color','r', 'FontWeight','bold', ...
         'HorizontalAlignment','center');
end
hold off;

