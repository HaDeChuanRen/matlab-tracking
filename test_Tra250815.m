clc; clear; close all;
rng(1);
mc_times = 10;
lw = 1;

num_frame = 50;
x_edge = 60e3;
y_edge = 40e3;
m_cbf = 361;
theta_vec = zeros(m_cbf, 1);
theta_vec(1 : floor(m_cbf / 2) + 1) = ...
    acos(linspace(1, -1, floor(m_cbf / 2) + 1))';
theta_vec(round(m_cbf / 2) : end) = ...
    acos(linspace(1, -1, floor(m_cbf / 2) + 1))' + pi;


n_range = 40e3;
k_targets = 100;
T_frame = n_range / 750;
vel_mean = 7;

lamb_clu = 200;
Pd_set = 1;
range_vec = 0 : n_range;

pos_ini = [-x_edge/2 + x_edge * rand(k_targets, 1), ...
    -y_edge/2 + y_edge * rand(k_targets, 1)];
vel_ini = [vel_mean * randn(k_targets, 1), vel_mean * randn(k_targets, 1)];
tar_track = zeros(k_targets, 2, num_frame);
mea_col = {};
x_cluedge = 3 * vel_mean * num_frame * T_frame + x_edge;
y_cluedge = 3 * vel_mean * num_frame * T_frame + y_edge;

% initialization
inxsig = 400;
inysig = 400;
Sigma_vec = [inxsig; inysig];
Nmax_hpy = 1;
N_mote = 4;
SPAMHT_plat = SPAMHTPlat(Nmax_hpy, Sigma_vec);

% HFM_last = [];
% timeStamp_last = 0;
% No_tol = 1;
% ini_Phyp = 1;
% time_cnt = 0;
% Tar_set = [];
% Clu_set = [];
% tar_his = [];
% K_show = [];
% Mea_idx = 1;
% Tra_idx = 1;
% P_D = 0.9;
% P_G = 1e-5;
% L_ite = 5;
% P_del = 0.1;
% 
% Tra_hyp = TraHypothesis(Tar_set, Clu_set, No_tol, ini_Phyp, time_cnt, ...
%     tar_his, Sigma_vec, P_G, L_ite, T_frame, Nmax_hpy, N_mote);
% Hypset_last = Tra_hyp;
% Max_hyp = 1;

tic;
for n_idx = 1 : num_frame
    % generate the position of targets
    tar_pos = pos_ini + (n_idx - 1) * T_frame * vel_ini;
    tar_track(:, :, n_idx) = tar_pos;
    tar_angle = wrapTo2Pi(atan2(tar_pos(:, 2), tar_pos(:, 1)));
    tar_mcbf = discretize(tar_angle, theta_vec);
    tar_det = [tar_mcbf, vecnorm(tar_pos, 2, 2)];

    num_clu = poissrnd(lamb_clu);
    clu_pos = [-x_cluedge/2 + x_cluedge * rand(num_clu, 1), ...
        -y_cluedge/2 + y_cluedge * rand(num_clu, 1)];
    clu_angle = wrapTo2Pi(atan2(clu_pos(:, 2), clu_pos(:, 1)));
    clu_mcbf = discretize(clu_angle, theta_vec);
    clu_det = [clu_mcbf, vecnorm(clu_pos, 2, 2)];

    mea_det = [tar_det; clu_det];

    % plot the detection result
%     h_fig1 = figure(1); clf(h_fig1);
%     plot(tar_det(:, 1), tar_det(:, 2), 'm*');
%     hold on; xlim([1 m_cbf]); ylim([0 n_range]);
%     plot(clu_det(:, 1), clu_det(:, 2), 'r*');
    
    angle_det = theta_vec(mea_det(:, 1));
    mea_pos = [mea_det(:, 2) .* cos(angle_det), ...
        mea_det(:, 2) .* sin(angle_det)];
    mea_col{end + 1} = mea_pos; %#ok<SAGROW>

    % data association and update
    [SPAMHT_plat, traj_all, target_all] = SPAMHT_plat.track(mea_pos);
%     Hypset_snap = [];
%     for h_idx = 1 : length(Hypset_last)
%         [Tra_hyp, All_hyp] = Hypset_last(h_idx).SPA_track(mea_pos);
%         Hypset_snap = [Hypset_snap; All_hyp]; %#ok<AGROW> 
%     end
%     
%     if length(Hypset_snap) > Max_hyp
%         Pro_hyp = zeros(length(Hypset_snap), 1);
%         for hp_idx = 1 : length(Hypset_snap)
%             Pro_hyp(hp_idx) = Hypset_snap(hp_idx).P_hyp;
%         end
%         [~, Psort_idx] = sort(Pro_hyp, 'descend');
%         Pro_prn = Pro_hyp(Psort_idx(1 : Max_hyp));
%         Pro_prn = Pro_prn / sum(Pro_prn);
%         Hypset_last = Hypset_snap(Psort_idx(1 : Max_hyp));
%         for h_idx = 1 : length(Hypset_last)
%             Hypset_last(h_idx).P_hyp = Pro_prn(h_idx);
%         end
%     else
%         Hypset_last = Hypset_snap;
%     end

    1;
    % pause(0.5);
end
toc;

h_fig2 = figure(2); hold on; grid on;
for k_idx = 1 : k_targets
    plot(squeeze(tar_track(k_idx, 1, :)), squeeze(tar_track(k_idx, 2, :)), ...
        'k--', 'LineWidth', lw);
end
for n_idx = 1 : num_frame
    plot(mea_col{n_idx}(:, 1), mea_col{n_idx}(:, 2), 'r*');
end

LineColors = hsv(10);
% Tra_hyp = Hypset_last(1);
% Tar_show = [Tra_hyp.Tar_set; Tra_hyp.Tra_his];
for t_idx = 1 : length(traj_all)
    plot(traj_all{t_idx}(:, 1), traj_all{t_idx}(:, 2), 'LineWidth', 1, ...
        'Color', LineColors(mod(t_idx, 10) + 1, :));
end

