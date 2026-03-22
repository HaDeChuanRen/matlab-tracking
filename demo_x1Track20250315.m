clear; clc; close all;
lw = 1;
fsz = 12;
msz = 5;
clim_up = 10;

addpath(genpath("E:\xmh_2025\2408syh\Matlab_syh"));

pathname = 'E:\xmh_2025\2408syh\Matlab_syh\Tracking\x1\testdata\output\';
filenames = dir([pathname '*.mat']);

% tracking parameter
N_trgu = [1000; 100];
CFAR_win = [ones(N_trgu(1) - N_trgu(2), 1); zeros(2 * N_trgu(2) + 1, 1);...
    ones(N_trgu(1) - N_trgu(2), 1)] / (2 * (N_trgu(1) - N_trgu(2)));
alpha_dB = 11;
cut_range = 100;
cut_cbf = 3;

t_vec = linspace(0, 2 * pi, 360);
Tmat = [1, 1; 0, 1];
inRsig = 100;
inAsig = 2;
P_D = 0.9;
P_G = 1e-5;
L_ite = 5;
P_del = 0.1;
T_ver = 52e3 / 750;


% initialization
HFM_last = [];
timeStamp_last = 0;
No_tol = 1;
ini_Phyp = 1;
time_cnt = 0;
Tar_set = [];
Clu_set = [];
tar_his = [];
K_show = [];
Mea_idx = 1;
Tra_idx = 1;
Sigma_vec = [inRsig; inAsig];
Nmax_hpy = 6;
N_mote = 4;

Tra_hyp = TraHypothesis(Tar_set, Clu_set, No_tol, ini_Phyp, time_cnt, ...
    tar_his, Sigma_vec, P_G, L_ite, T_ver, Nmax_hpy, N_mote);
Hypset_last = Tra_hyp;
Max_hyp = 1;

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
    % mw_hfm = mw_data_hfm.';
    [N_range, N_cbf] = size(mw_hfm);
    Rist_min = 1;
    Rist_max = N_range - N_trgu(1);
    A_vec = acosd(linspace(-1, 1, N_cbf / 2));
%     [R_grids, Cbf_grids] = meshgrid([- A_vec, A_vec], ...
%         (Rist_min : Rist_max) / 1000);
    [R_grids, Cbf_grids] = meshgrid(1 : N_cbf, (Rist_min : Rist_max) / 1000);

    [RA_grids, A_grids] = meshgrid([- A_vec, A_vec], (1 : N_range) / 1000);
    h_fig3 = figure(3);
    pcolor(R_grids, Cbf_grids, mw_hfm(Rist_min : Rist_max, :));
    colormap(jet); shading flat; clim([0 clim_up]);
    xlabel('波束');
    ylabel('距离/km');

    h_fig1 = figure(1); h_pclr1 = pcolor(R_grids, Cbf_grids, ...
        mw_hfm(Rist_min : Rist_max, :));
    set(h_pclr1, 'EdgeColor', 'none'); colormap(jet);
    clim([0 clim_up]);
    hold on; fig_title = ['period = ', num2str(ifile)];
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
        det_nc = [r_det, nc_idx * ones(length(r_det), 1)];
        mea_det = [mea_det; det_nc]; %#ok<*AGROW>
        mea_snr = [mea_snr; mw_snr(r_det)];
    end

    % plot the result of detection
    h_fig2 = figure(2); clf(h_fig2); 
    h_pclr2 = pcolor(R_grids, Cbf_grids, mw_hfm(Rist_min : Rist_max, :));
    set(h_pclr2, 'EdgeColor', 'none'); colormap(jet);
    clim([0 clim_up]);
    hold on; title(fig_title);
    plot(mea_det(:, 2), mea_det(:, 1) / 1000, 'w*', 'Markersize', msz);

    mea_temp = zeros(size(mea_det, 1), 2);
    mea_temp(:, 2) = mea_det(:, 2) * 30;
    mea_temp(:, 1) = mea_det(:, 1);
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
        size_r = max(mea_det(idx_clun, 1)) - min(mea_det(idx_clun, 1));
        size_cbf = max(mea_det(idx_clun, 2)) - min(mea_det(idx_clun, 2));
        % if size_r / size_cbf < 1 / 3 % 用强度和横向尺寸的比例来计算
        %     idx_del = [idx_del; idx_clun];
        %     continue;
        % end
        N_clu = N_clu + 1;
        clu_n = [clu_n; n_idx];
        mea_clun = mean(mea_det(idx_clun, :), 1);
        ellipse_pos = mea_clun' + [(size_r / 2 + 20) * sin(t_vec); ...
            (size_cbf / 2 + 0.5) * cos(t_vec)];
        figure(h_fig2); plot(ellipse_pos(2, :), ellipse_pos(1, :) / 1000,  ...
            'r-','LineWidth', lw);
        mea_clu(end + 1, :) = mea_clun;
        meares_cell = {string(timeStamp_last), mea_clun(1), mea_clun(2)};
        mea_w = ['A', num2str(Mea_idx)];
        writecell(meares_cell, 'detection_res.xlsx', 'Sheet', 'Sheet1', ...
            'Range', mea_w);
        Mea_idx = Mea_idx + 1;
    end
    figure(h_fig2); hold off;
    mea_det(idx_del, :) = []; %#ok<*SAGROW> 
    res_clu(idx_del) = [];
    mea_snr(idx_del) = [];

    % data association and update
    Hypset_snap = [];
    for h_idx = 1 : length(Hypset_last)
        [Tra_hyp, All_hyp] = Hypset_last(h_idx).SPA_track(mea_clu);
        Hypset_snap = [Hypset_snap; All_hyp];
    end
    
    
    if length(Hypset_snap) > Max_hyp
        Pro_hyp = zeros(length(Hypset_snap), 1);
        for hp_idx = 1 : length(Hypset_snap)
            Pro_hyp(hp_idx) = Hypset_snap(hp_idx).P_hyp;
        end
        [~, Psort_idx] = sort(Pro_hyp, 'descend');
        Pro_prn = Pro_hyp(Psort_idx(1 : Max_hyp));
        Pro_prn = Pro_prn / sum(Pro_prn);
        Hypset_last = Hypset_snap(Psort_idx(1 : Max_hyp));
        for h_idx = 1 : length(Hypset_last)
            Hypset_last(h_idx).P_hyp = Pro_prn(h_idx);
        end
    else
        Hypset_last = Hypset_snap;
    end

    Tar_set = Hypset_last(1).Tar_set;
    K_tar = length(Tar_set);
    K_show = [];
    for k_idx = 1 : K_tar
        if Tar_set(k_idx).P_s >= 0.6
            K_show(end + 1) = k_idx;
        end
    end

    % plot the tracking result
    Num_show = length(K_show);
    det_show = zeros(Num_show, 2);
    % Grade_show = zeros(Num_show, 1);
    for det_idx = 1 : Num_show
        k_idx = K_show(det_idx);
        r_k = Tar_set(k_idx).State_vec(1) / 1000;
        an_k = Tar_set(k_idx).State_vec(2);
        det_show(det_idx, :) = [an_k, r_k];
        % Grade_show(det_idx) = Tar_set(k_idx).P_s * Tar_set(k_idx).Num_vala;
        figure(h_fig1); plot(an_k, r_k, 'wo', 'MarkerSize', msz);
        rhis_k = [Tar_set(k_idx).his_state(:, 1) / 1000; r_k];
        ahis_k = [Tar_set(k_idx).his_state(:, 2); an_k];
        figure(h_fig1); plot(ahis_k, rhis_k, 'w-', 'Linewidth', lw);
        str_tar = {num2str(Tar_set(k_idx).No_tar), ...
            [num2str(sum(Tar_set(k_idx).his_Val)), ',', ...
            num2str(Tar_set(k_idx).P_s)]};
        text(an_k + 3, r_k - 0.3, str_tar, 'Color','white','FontSize', 6);

        trares_cell = {string(timeStamp_last), Tar_set(k_idx).No_tar, ...
            string(Tar_set(k_idx).his_t(1)), r_k, an_k};
        range_w = ['A', num2str(Tra_idx)];
        writecell(trares_cell, 'Tracking_res.xlsx', 'Sheet', 'Sheet1', ...
            'Range', range_w);
        Tra_idx = Tra_idx + 1;
    end

    h_fig4 = figure(4); h_pclr4 = pcolor(R_grids, Cbf_grids, ...
        mw_hfm(Rist_min : Rist_max, :));
    set(h_pclr4, 'EdgeColor', 'none'); colormap(jet);
    clim([0 clim_up]);
    hold on; fig_title = ['period = ', num2str(ifile)];
    title(fig_title);
    for det_idx = 1 : Num_show
        k_idx = K_show(det_idx);
        r_k = Tar_set(k_idx).State_vec(1) / 1000;
        an_k = Tar_set(k_idx).State_vec(2);
        vr_k = Tar_set(k_idx).Vel_vec(1) / T_ver;

        if (abs(vr_k) > 0.5)
            figure(h_fig4); plot(an_k, r_k, 'r*', 'MarkerSize', msz);
            str_tar = {num2str(Tar_set(k_idx).No_tar), num2str(vr_k, '%.2f')};
            text(an_k + 3, r_k - 0.3, str_tar, 'Color','white','FontSize', 6);
        else
            figure(h_fig4); plot(an_k, r_k, 'g*', 'MarkerSize', msz);
            str_tar = {num2str(Tar_set(k_idx).No_tar), num2str(vr_k, '%.2f')};
            text(an_k + 3, r_k - 0.3, str_tar, 'Color','white','FontSize', 6);
        end
    end

    % frame_fig1 = getframe(h_fig1);
    % writeVideo(writer_fig1, frame_fig1);
    h_fig1; hold off;
    h_fig4; hold off;
    if ifile == 9
        1;
    end
end
% close(writer_fig1);