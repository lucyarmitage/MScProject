clear;

ref_mat   = 'Prostate_Ref_T1_T2.mat';
roi_file  = 'Prostate_Reference.roi';   % ImageJ ROI

T1_RANGE = [0 1500];
T2_RANGE = [0 250];

% Load reference T1 and T2 maps
ref = load(ref_mat);
T1_ref = double(ref.T1map);
T2_ref = double(ref.T2map);
[H, W] = size(T1_ref);


% Read ImageJ ROI
roi_struct = ReadImageJROI(roi_file);
% If multiple ROIs returned, pick the first 'Oval'
if numel(roi_struct) > 1
    idx_oval = find(arrayfun(@(r) isfield(r,'strType') && strcmpi(r.strType,'Oval'), roi_struct), 1);
    if isempty(idx_oval), idx_oval = 1; end
    roi_struct = roi_struct(idx_oval);
end

% Build ellipse mask on reference image grid
ellipseROI.strType = 'Oval';
ellipseROI.vnRectBounds = roi_struct.vnRectBounds;
mask = ellipse2mask(ellipseROI, [H, W]);   % logical(H,W)


% ROI on top of reference T1 & T2
figure('Name','T1 reference + ROI');
show_overlay(T1_ref, mask, T1_RANGE, 'T1 reference with ROI');

figure('Name','T2 reference + ROI');
show_overlay(T2_ref, mask, T2_RANGE, 'T2 reference with ROI');

% Compute ROI means
T1_vals_ref   = T1_ref(mask);
T2_vals_ref   = T2_ref(mask);
mean_T1_ref   = mean(T1_vals_ref,   'omitnan');
mean_T2_ref   = mean(T2_vals_ref,   'omitnan');
fprintf('Mean T1 (reference):      %.2f ms\n', mean_T1_ref);
fprintf('Mean T2 (reference):      %.2f ms\n', mean_T2_ref);


function show_overlay(img2d, mask, clim, titleStr)

    img2d = double(img2d);
    mask  = logical(mask);

    % Base image
    imagesc(img2d, clim);
    axis image off;
    set(gca,'YDir','normal');      % standard orientation
    colorbar;
    title(titleStr);

    hold on;

    cyanRGB = cat(3, zeros(size(mask)), ones(size(mask)), ones(size(mask))); % [R G B]
    hFill = image(cyanRGB);
    set(hFill, 'AlphaData', 0.25 * mask, 'HitTest','off');
    B = bwperim(mask);
    [y, x] = find(B);
    plot(x, y, 'c', 'LineWidth', 1.5);

    hold off;
    drawnow;
end

function mask = ellipse2mask(roi, imageSize)
% ellipse2mask: Creates a binary mask for an ImageJ ellipse ROI
% Usage:
%   mask = ellipse2mask(roi, [height, width]);
%   mask = ellipse2mask(roi, imageStruct);  % where imageStruct.nHeight, nWidth

    if isnumeric(imageSize) && numel(imageSize) == 2
        imageStruct.nHeight = imageSize(1);
        imageStruct.nWidth  = imageSize(2);
        imageSize = imageStruct;
    end

    if ~isfield(roi, 'strType') || ~strcmpi(roi.strType, 'Oval')
        error('ellipse2mask: Input ROI must be of type ''Oval''.');
    end
    if ~isfield(roi, 'vnRectBounds') || numel(roi.vnRectBounds) ~= 4
        error('ellipse2mask: ROI must contain vnRectBounds with 4 elements.');
    end

    % Extract bounding box
    x1 = roi.vnRectBounds(1);
    y1 = roi.vnRectBounds(2);
    x2 = roi.vnRectBounds(3);
    y2 = roi.vnRectBounds(4);

    % Create meshgrid
    [xx, yy] = meshgrid(1:imageSize.nWidth, 1:imageSize.nHeight);

    % Compute center and radii
    xc = (x1 + x2) / 2;
    yc = (y1 + y2) / 2;
    rx = max((x2 - x1) / 2, eps);
    ry = max((y2 - y1) / 2, eps);

    % Ellipse equation (x-xc)^2/rx^2 + (y-yc)^2/ry^2 <= 1
    mask = ((xx - xc).^2) / rx^2 + ((yy - yc).^2) / ry^2 <= 1;
    mask = logical(mask);
end
