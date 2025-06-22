folder = 'C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/';
load(fullfile(folder, 'recon_results101.mat')); 

% Show the image 
magImg = abs(svd_images(:,:,1));
magImg = magImg / max(magImg(:));

f = figure; imshow(magImg, []); title('Draw the first circle to fix radius');

% Draw first circle
h = drawcircle('Color','r');
wait(h);  % Waits to finish adjusting
baseRadius = h.Radius;

nVials = 14;
circles = gobjects(nVials,1);
centres = zeros(nVials,2);
radii   = baseRadius * ones(nVials,1);
circles(1) = h;  % Store first circle

% Click centre positions for remaining vials
title(sprintf('Click centre of remaining %d vials', nVials - 1));
[x, y] = ginput(nVials - 1);

for k = 2:nVials
    circles(k) = drawcircle('Center', [x(k-1), y(k-1)], ...
                            'Radius', baseRadius, 'Color','r');
end

% Create "Done" button
btn = uicontrol('Style', 'pushbutton', 'String', 'Done', ...
                'Position', [20 20 80 30], ...
                'Callback', 'uiresume(gcbf);');

title('Adjust circles if needed, then click "Done"');

% Wait until Done is pressed
uiwait(f); 

% Extract centres and radii (safe!)
for k = 1:nVials
    centres(k,:) = circles(k).Center;
    radii(k)     = circles(k).Radius;
end

% Close figure and save result
close(f);
save('vial_mask_geometry.mat', 'centres', 'radii');
disp('Saved vial positions to vial_mask_geometry.mat');
