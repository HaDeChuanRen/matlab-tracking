clear; clc; close all;

% plot parameters
lw = 1;
fsz = 12;
msz = 5;
clim_up = 10;
t_vec = linspace(0, 2 * pi, 360);

% file path
pathname = 'E:\xmh_2025\2408syh\Matlab_syh\Tracking\x1\testdata\output\';
filenames = dir([pathname '*.mat']);

% detection initialization
N_tra = [500, 0];
N_gua = [50, 0];
alpha_dB = 11;
CA_det = CACFARdet(N_tra, N_gua);

% tracking initialization
inxsig = 100;
inysig = 2;
Sigma_vec = [inxsig; inysig];
Nmax_hpy = 1;
% N_mote = 4;
SPAMHT_plat = SPAMHTPlat(Nmax_hpy, Sigma_vec);
HFM_last = [];

% Tmat = [1, 1; 0, 1];
% inRsig = 100;
% inAsig = 2;
% P_D = 0.9;
% P_G = 1e-5;
% L_ite = 5;
% P_del = 0.1;
% T_ver = 52e3 / 750;

% initialization
% % timeStamp_last = 0;
% No_tol = 1;
% ini_Phyp = 1;
% time_cnt = 0;
% Tar_set = [];
% Clu_set = [];
% tar_his = [];
% K_show = [];
% Mea_idx = 1;
% Tra_idx = 1;
% Sigma_vec = [inRsig; inAsig];
% Nmax_hpy = 6;
% N_mote = 4;
% 
% Tra_hyp = TraHypothesis(Tar_set, Clu_set, No_tol, ini_Phyp, time_cnt, ...
%     tar_his, Sigma_vec, P_G, L_ite, T_ver, Nmax_hpy, N_mote);
% Hypset_last = Tra_hyp;
% Max_hyp = 1;

for ifile = 1:length(filenames)
    %% obtain and compensate HFM
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

    %% plot HFM
    mw_hfm = (HFM_move') .^ 2;
    [N_range, N_cbf] = size(mw_hfm);
    Rist_min = 1e4;
    Rist_max = N_range - N_tra(1);
    A_vec = acosd(linspace(-1, 1, N_cbf / 2));

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


    %% Detection and Clustering
    mea_det = CA_det.detection2D(mw_hfm(Rist_min : Rist_max, :), alpha_dB);
    mea_det(:, 1) = mea_det(:, 1) + Rist_min - 1;

    % plot the result of detection
    h_fig2 = figure(2); clf(h_fig2); 
    h_pclr2 = pcolor(R_grids, Cbf_grids, mw_hfm(Rist_min : Rist_max, :));
    set(h_pclr2, 'EdgeColor', 'none'); colormap(jet);
    clim([0 clim_up]);
    hold on; title(fig_title);
    plot(mea_det(:, 2), mea_det(:, 1) / 1000, 'w*', 'Markersize', msz);
    [mea_clu, size_clu, res_clu] = cluster_link(mea_det, [1 30]);
    N_clu = size(mea_clu, 1);

    for n_idx = 1 : N_clu
        mea_clun = mea_clu(n_idx, :);
        ellipse_pos = mea_clun' + [(size_clu(n_idx, 1) / 2 + 20) * sin(t_vec); ...
            (size_clu(n_idx, 2) / 2 + 0.5) * cos(t_vec)];
        figure(h_fig2); plot(ellipse_pos(2, :), ellipse_pos(1, :) / 1000,  ...
            'r-','LineWidth', lw);
    end
    figure(h_fig2); hold off;

    %% data association and tracking
    [SPAMHT_plat, ~, Tar_set] = SPAMHT_plat.track(mea_clu);

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
    for det_idx = 1 : Num_show
        k_idx = K_show(det_idx);
        r_k = Tar_set(k_idx).State_vec(1) / 1000;
        an_k = Tar_set(k_idx).State_vec(2);
        det_show(det_idx, :) = [an_k, r_k];
        figure(h_fig1); plot(an_k, r_k, 'wo', 'MarkerSize', msz);
        rhis_k = [Tar_set(k_idx).his_state(:, 1) / 1000; r_k];
        ahis_k = [Tar_set(k_idx).his_state(:, 2); an_k];
        figure(h_fig1); plot(ahis_k, rhis_k, 'w-', 'Linewidth', lw);
        str_tar = {num2str(Tar_set(k_idx).No_tar), ...
            [num2str(sum(Tar_set(k_idx).his_Val)), ',', ...
            num2str(Tar_set(k_idx).P_s)]};
        text(an_k + 3, r_k - 0.3, str_tar, 'Color','white','FontSize', 6);
    end

    h_fig1; hold off;

    if ifile == 9
        1;
    end

end
% close(writer_fig1);