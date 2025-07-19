clear;

recon_folder = 'C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/recon_results';
load(fullfile(recon_folder, '10mm_501_short.mat'));  % Loads svd_images, D, dict0, idx, R, match_images, idx2, c, N

gs_data = load('T1T2_results.mat');
t1_gs = gs_data.T1map;
t2_gs = gs_data.T2map;

x = 61;   
y = 50;  
x_flipped = N - x + 1;  % reverse fliplr
y_flipped = N - y + 1;  % reverse flip
voxel_idx = sub2ind([N N], y_flipped, x_flipped);


voxel_signal = match_images(voxel_idx, :);
sim_vector = abs(voxel_signal * conj(D.magnetization));  


[sorted_vals, sorted_idx] = sort(sim_vector, 'descend');
top5_idx = sorted_idx(1:5);
top5_sim = sorted_vals(1:5);
top5_T1 = D.lookup_table(top5_idx, 1);
top5_T2 = D.lookup_table(top5_idx, 2);


if ~isequal(size(t1_gs), [N N])
    t1_gs = imresize(t1_gs, [N N]);
    t2_gs = imresize(t2_gs, [N N]);
end

gs_T1 = t1_gs(y, x);
gs_T2 = t2_gs(y, x);

fprintf('Top 5 dictionary matches at pixel (x=%d, y=%d):\n', x, y);
fprintf('GS T1: %.1f ms, GS T2: %.1f ms\n', gs_T1, gs_T2);
fprintf('%-5s %-10s %-10s %-10s\n', '#', 'T1 (ms)', 'T2 (ms)', 'Similarity');
for i = 1:5
    fprintf('%-5d %-10.1f %-10.1f %-10.4f\n', ...
        i, top5_T1(i), top5_T2(i), top5_sim(i));
end
