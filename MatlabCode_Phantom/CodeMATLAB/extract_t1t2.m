clear;

folder = 'C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/recon_results';
addpath(genpath('C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/CodeMATLAB/utils'));
load(fullfile(folder, 'new_10mm_501_short.mat'));
load('vial_mask_geometry.mat');  % centres, radii

Qmaps = reshape(D.lookup_table(idx2,:),[[N N 1], size(D.lookup_table,2)]);
Qmaps = cat(numel(size(Qmaps)),Qmaps,reshape(c ./ D.normalization(idx2).',[N N]));
% PD    = c ./ D.normalization(idx2).';
% PD    = reshape(PD, [N N]);

% Do some background masking
mask = sum(abs(reshape(svd_images,[[N N] R])),3);
mask(mask < .4* max(mask(:))) = 0;
mask (mask ~= 0) = 1;
Qmaps = bsxfun(@times,Qmaps,mask);

Qmaps = flip (Qmaps);
% figure; imshow3(fliplr(Qmaps(:,:,1,1)), [0 max(max(Qmaps(:,:,1,1)))])
% figure; imshow3(fliplr(Qmaps(:,:,1,2)), [0 max(max(Qmaps(:,:,1,2)))])
% max(max(Qmaps(:,:,1,1)))
% max(max(Qmaps(:,:,1,2)))

figure; imshow3(fliplr(Qmaps(:,:,1,1)), [0 3000]);
figure; imshow3(fliplr(Qmaps(:,:,1,2)), [0 1500]); 
% figure; imshow3(fliplr(Qmaps(:,:,1,1)), [0 3000]); colormap hot;
% figure; imshow3(fliplr(Qmaps(:,:,1,2)), [0 1500]); colormap turbo;


magImg = fliplr(Qmaps(:,:,1,1));
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

figure; 
hImg = imagesc(magImg); axis image; colormap hot; colorbar;
title('Click on pixel to view T1/T2 values');

dcm = datacursormode(gcf);
set(dcm, 'Enable', 'on', ...
         'DisplayStyle', 'datatip', ...
         'SnapToDataVertex', 'on', ...
         'UpdateFcn', {@myTooltipFunction, T1, T2});


disp('Click');


function txt = myTooltipFunction(~, event_obj, T1, T2)
    pos = get(event_obj, 'Position');
    x = round(pos(1));
    y = round(pos(2));

    if y >= 1 && y <= size(T1, 1) && x >= 1 && x <= size(T1, 2)
        t1_val = T1(y, x);
        t2_val = T2(y, x);
        txt = {
            ['X: ', num2str(x)], ...
            ['Y: ', num2str(y)], ...
            ['T1: ', num2str(t1_val, '%.1f'), ' ms'], ...
            ['T2: ', num2str(t2_val, '%.1f'), ' ms']
        };
    else
        txt = {'Out of bounds'};
    end
end















% --- Pick a pixel and report its index ---
[x_click, y_click] = ginput(1);
x = round(x_click);
y = round(y_click);
pixel_idx = sub2ind([N, N], y, x);  % Make sure N matches actual image size

fprintf('You selected pixel (x = %d, y = %d)\n', x, y);
fprintf('Linear index (pixel_idx) = %d\n', pixel_idx);

