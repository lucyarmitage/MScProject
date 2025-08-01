function interactive_voxel_viewer_callback()
    close all;
    clear;

    folder = 'C:/Users/lucya/MSC_PROJECT/MatlabCode_Phantom/recon_results';
    dict_files = {
        '10mm_101_short.mat'
        '10mm_101_short_b1.mat'
    };
    dict_labels = {
        'normal'
        'b1'
    };

    data0 = load(fullfile(folder, dict_files{1}));
    svd_images = data0.svd_images;
    Qmaps = data0.Qmaps;
    N = size(svd_images, 1);
    
    % Gold standard 
    gs_data = load('T1T2_results.mat');
    t1_gs = imresize(gs_data.T1map, [N N]);
    t2_gs = imresize(gs_data.T2map, [N N]);

    map_to_show = Qmaps(:,:,1,1);
    fig = figure(1); clf;
    imshow(fliplr(map_to_show), [0 3000]);

    set(fig, 'WindowButtonDownFcn', ...
        @(src, event) on_voxel_click(src, event, folder, dict_files, dict_labels, N, Qmaps, svd_images, t1_gs, t2_gs));
end

function on_voxel_click(src, ~, folder, dict_files, dict_labels, N, Qmaps, svd_images, t1_gs, t2_gs)
    coords = get(gca, 'CurrentPoint');
    x_click = round(coords(1,1));
    y_click = round(coords(1,2));

    x = N - x_click + 1;
    y = N - y_click + 1;
    pixel_idx = sub2ind([N, N], y, x);

    % dictionary
    base = load(fullfile(folder, dict_files{1}));
    D = base.D;
    idx = base.idx;
    idx2 = base.idx2;
    dict0 = base.dict0;

    % measured signal
    compressed = squeeze(svd_images(y, x, :));
    meas_sig = D.u * compressed;

    % dictionary matched result
    matched_idx = idx2(pixel_idx);
    T1_match = idx(matched_idx, 1);
    T2_match = idx(matched_idx, 2);

    % gold standard
    gs_T1 = t1_gs(y_click, x_click);
    gs_T2 = t2_gs(y_click, x_click);

    % Find closest dictionary entry to ground truth
    dists = sqrt((idx(:,1) - gs_T1).^2 + (idx(:,2) - gs_T2).^2);
    [~, gs_idx] = min(dists);
    gs_T1_dict = idx(gs_idx, 1);
    gs_T2_dict = idx(gs_idx, 2);

    fprintf('\n Voxel (%d,%d) \n', x_click, y_click);
    fprintf('GS T1 = %.1f, GS T2 = %.1f\n', gs_T1, gs_T2);
    fprintf('Matched T1 = %.1f, T2 = %.1f\n', T1_match, T2_match);
    fprintf('Closest dict entry to GS: T1 = %.1f, T2 = %.1f\n', gs_T1_dict, gs_T2_dict);

    %% Plot

    figure(2); clf; f_abs = gcf;
    figure(3); clf; f_real = gcf;
    figure(4); clf; f_imag = gcf;

    legend_abs = {};
    legend_real = {};
    legend_imag = {};
    colors = lines(numel(dict_files));

    for i = 1:numel(dict_files)
        data = load(fullfile(folder, dict_files{i}));
        dict0 = data.dict0;

        dict_sig = conj(dict0(:, matched_idx));
        dict_sig_scaled = dict_sig * (norm(meas_sig) / norm(dict_sig));
        dict_sig_aligned = dict_sig_scaled * exp(-1i * angle(sum(conj(meas_sig) .* dict_sig_scaled)));

        gs_dict_sig = dict0(:, gs_idx);
        gs_dict_sig_scaled = gs_dict_sig * (norm(meas_sig) / norm(gs_dict_sig));
        gs_dict_sig_aligned = gs_dict_sig_scaled * exp(-1i * angle(sum(conj(meas_sig) .* gs_dict_sig_scaled)));

        figure(f_abs);
        plot(abs(dict_sig_aligned), '--', 'Color', colors(i,:), 'LineWidth', 1.2); hold on;
        plot(abs(gs_dict_sig_aligned), '-', 'Color', colors(i,:), 'LineWidth', 1.2);
        legend_abs{end+1} = sprintf('%s Matched', dict_labels{i});
        legend_abs{end+1} = sprintf('%s GS dict', dict_labels{i});

        figure(f_real);
        plot(real(dict_sig_aligned), '--', 'Color', colors(i,:), 'LineWidth', 1.2); hold on;
        plot(real(gs_dict_sig_aligned), '-', 'Color', colors(i,:), 'LineWidth', 1.2);
        legend_real{end+1} = sprintf('%s Matched', dict_labels{i});
        legend_real{end+1} = sprintf('%s GS dict', dict_labels{i});

        figure(f_imag);
        plot(imag(dict_sig_aligned), '--', 'Color', colors(i,:), 'LineWidth', 1.2); hold on;
        plot(imag(gs_dict_sig_aligned), '-', 'Color', colors(i,:), 'LineWidth', 1.2);
        legend_imag{end+1} = sprintf('%s Matched', dict_labels{i});
        legend_imag{end+1} = sprintf('%s GS dict', dict_labels{i});
    end

    % measured signal
    figure(f_abs);
    plot(abs(meas_sig), 'k', 'LineWidth', 2);
    title(sprintf('Abs: Voxel (%d,%d)', x_click, y_click));
    xlabel('Time'); ylabel('|Signal|'); grid on;
    legend_abs{end+1} = 'Measured';
    legend(legend_abs, 'Location', 'eastoutside');

    figure(f_real);
    plot(real(meas_sig), 'k', 'LineWidth', 2);
    title(sprintf('Real: Voxel (%d,%d)', x_click, y_click));
    xlabel('Time'); ylabel('Real'); grid on;
    legend_real{end+1} = 'Measured';
    legend(legend_real, 'Location', 'eastoutside');

    figure(f_imag);
    plot(imag(meas_sig), 'k', 'LineWidth', 2);
    title(sprintf('Imag: Voxel (%d,%d)', x_click, y_click));
    xlabel('Time'); ylabel('Imag'); grid on;
    legend_imag{end+1} = 'Measured';
    legend(legend_imag, 'Location', 'eastoutside');
end
