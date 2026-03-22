clear; clc; close all;
pathname = 'D:\xmh_2025\202408syh\Matlab_syh\Tracking\x1\testdata\output\';
filenames = dir([pathname '*.mat']);
lw = 1;
fsz = 12;
msz = 5;

% tracking parameter
N_trgu = [1000; 100];
CFAR_win = [ones(N_trgu(1) - N_trgu(2), 1); zeros(2 * N_trgu(2) + 1, 1);...
    ones(N_trgu(1) - N_trgu(2), 1)] / (2 * (N_trgu(1) - N_trgu(2)));
alpha_dB = 11;
cut_range = 100;
cut_cbf = 3;
T_ver = 120;
Tar_set = [];

t_vec = linspace(0, 2 * pi, 360);
Tmat = [1, 1; 0, 1];
inRsig = 200;
inAsig = 3;
P_D = 0.9;
P_G = 1e-6;
L_ite = 5;

% create mp4 file
% mp4_fig1 = 'mp4_save\track_x1.mp4';
% writer_fig1 = VideoWriter(mp4_fig1, 'MPEG-4');
% writer_fig1.FrameRate = 3;
% open(writer_fig1);


% initialization
HFM_last = [];
timeStamp_last = 0;
No_tol = 1;
row_w = 2;

for ifile = 1:length(filenames)
    load([pathname filenames(ifile).name]);

    if ~isempty(HFM_last)
        timeGap = seconds(command.timeStamp - timeStamp_last);
        timeStamp_last = datetime(command.timeStamp);
        clear("HFM_last");
        HFM_last = HFM;
        clear("HFM");
    else
        timeStamp_last = command.timeStamp;
        clear("HFM_last");
        HFM_last = HFM;
        clear("HFM");
        continue;
    end

    platSpeed = command.platSpeed / 10;
    platX = timeGap * platSpeed / 1.9438;
    cosSideAngle = 1 : -2/327 : -1;
    rX = platX.*cosSideAngle;
    HFM_move = zeros(size(HFM_last));
    for ibeam = 1:328
        HFM_move(ibeam,:) = circshift(HFM_last(ibeam,:),-floor(rX(ibeam)));
    end


    mw_hfm = (HFM_move') .^ 2;
    [N_range, N_cbf] = size(mw_hfm);
    Rist_min = 10000;
    Rist_max = N_range - N_trgu(1);
    [R_grids, Cbf_grids] = meshgrid(1 : N_cbf, Rist_min : Rist_max);

    A_vec = acosd(linspace(-1, 1, N_cbf / 2));
    [RA_grids, A_grids] = meshgrid([- A_vec, A_vec], (1 : N_range) / 1000);
    h_fig3 = figure(3); pcolor(RA_grids, A_grids, mw_hfm);
    colormap(jet); shading flat; clim([0 10]);
    xlabel('方位角/°');
    ylabel('距离/km');

    h_fig1 = figure(1); h_pclr1 = pcolor(R_grids, Cbf_grids, ...
        mw_hfm(Rist_min : Rist_max, :));
    set(h_pclr1, 'EdgeColor', 'none'); colormap(jet); clim([0 10]); hold on;
    fig_title = ['period = ', num2str(ifile)];
    title(fig_title);

    mea_det = [];
    mea_snr = [];
    for nc_idx = 1 : N_cbf
        mw_sigma = conv(CFAR_win, mw_hfm(:, nc_idx));
        mw_snr = 10 * log10(mw_hfm(:, nc_idx) ./ mw_sigma(N_trgu(1) + ...
            (1 : N_range)));
        r_det = find(mw_snr > alpha_dB);
        r_det(r_det < Rist_min) = [];
        r_det(r_det > Rist_max) = [];
        det_nc = [nc_idx * ones(length(r_det), 1), r_det];
        mea_det = [mea_det; det_nc]; %#ok<*AGROW>
        mea_snr = [mea_snr; mw_snr(r_det)];
    end

    % plot the result of detection
    h_fig2 = figure(2); clf(h_fig2); 
    h_pclr2 = pcolor(R_grids, Cbf_grids, mw_hfm(Rist_min : Rist_max, :));
    set(h_pclr2, 'EdgeColor', 'none'); colormap(jet); clim([0 10]); hold on;
    title(fig_title);
    plot(mea_det(:, 1), mea_det(:, 2), 'w*', 'Markersize', msz);

    mea_temp = zeros(size(mea_det, 1), 2);
    mea_temp(:, 1) = mea_det(:, 1) * 30;
    mea_temp(:, 2) = mea_det(:, 2);
    Dis_tar = pdist(mea_temp);
    Fre_link = linkage(Dis_tar, 'single'); %figure; dendrogram(Fre_link)
    res_clu = cluster(Fre_link, 'Cutoff', 100, 'Criterion', 'distance');
    N_max = max(res_clu);

    N_clu = 0;
    idx_del = [];
    clu_n = [];
    mea_clu = [];
    for n_idx = 1 : N_max
        idx_clun = find(res_clu == n_idx);
        if isempty(idx_clun)
            continue;
        end
        size_r = max(mea_det(idx_clun, 2)) - min(mea_det(idx_clun, 2));
        size_cbf = max(mea_det(idx_clun, 1)) - min(mea_det(idx_clun, 1));
        if size_r / size_cbf < 1 / 3 % 用强度和横向尺寸的比例来计算
            idx_del = [idx_del; idx_clun];
            continue;
        end
        N_clu = N_clu + 1;
        clu_n = [clu_n; n_idx];
        mea_clun = mean(mea_det(idx_clun, :), 1);
        ellipse_pos = mea_clun' + [(size_cbf / 2 + 0.5) * cos(t_vec); ...
            (size_r / 2 + 20) * sin(t_vec)];
        figure(h_fig2); plot(ellipse_pos(1, :), ellipse_pos(2, :),  ...
            'r-','LineWidth', lw);
        mea_clu(end + 1, :) = mea_clun;
    end
    figure(h_fig2); hold off;
    mea_det(idx_del, :) = []; %#ok<*SAGROW> 
    res_clu(idx_del) = [];
    mea_snr(idx_del) = [];

    % data association and update
    [Tar_set, Promat_a, km_ass, K_show] = SPA_Track(Tar_set, mea_clu, ...
        timeStamp_last, P_G, inAsig, inRsig);

    % plot the tracking result
    Num_show = length(K_show);
    det_show = zeros(Num_show, 2);
    Grade_show = zeros(Num_show, 1);
    for det_idx = 1 : Num_show
        k_idx = K_show(det_idx);
        r_k = Tar_set(k_idx).state_R(1);
        an_k = Tar_set(k_idx).state_A(1);
        det_show(det_idx, :) = [an_k, r_k];
        Grade_show(det_idx) = Tar_set(k_idx).P_s * Tar_set(k_idx).Num_vala;
        figure(h_fig1); plot(an_k, r_k, 'wo', 'MarkerSize', msz);
        rhis_k = [Tar_set(k_idx).Range_his, r_k];
        ahis_k = [Tar_set(k_idx).Angle_his, an_k];
        figure(h_fig1); plot(ahis_k, rhis_k, 'w-', 'Linewidth', lw);
        str_tar = {num2str(Tar_set(k_idx).No_tar), ...
            [num2str(Tar_set(k_idx).Num_vala), ',', ...
            num2str(Tar_set(k_idx).P_s)]};
        text(an_k + 3, r_k - 3, str_tar, 'Color','white','FontSize', 6);

        trares_cell = {string(timeStamp_last), Tar_set(k_idx).No_tar, ...
            string(Tar_set(k_idx).t_his(1)), r_k, an_k};
        range_w = ['A', num2str(row_w)];
        writecell(trares_cell, 'Tracking_res.xlsx', 'Sheet', 'Sheet1', ...
            'Range', range_w);
        row_w = row_w + 1;
    end

    frame_fig1 = getframe(h_fig1);
    % writeVideo(writer_fig1, frame_fig1);
    hold off;
end
% close(writer_fig1);