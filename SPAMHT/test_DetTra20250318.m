clc; clear; close all;

lw = 1.5;
msz = 8;
fsz = 12;

addpath(genpath('E:\xmh_2025\202502PasTracking\'));

% load('Data_2x6B_3_BF.mat');
% Tra_mat = cbf3;
% load('MultiBeam_20250226020525.mat');
% Tra_mat = gather(MB)';

load('FZ.mat');
Tra_mat = BTR2;

% load('40UTY_20250314_160512_SubRcbEck_100_200Hz_16s_misgrate0s.mat');
% Tra_mat=squeeze(result(:,1:2000,2)).';


% load('XuMenghuai_2_20250228_11_1.mat', 'Tra_mat');
% Tra_mat = res_cbf;

[Nsnap, NB] = size(Tra_mat);
figure; pcolor(Tra_mat); shading flat; colormap jet; % clim([0, 6000]);

% CFAR detector parameter
num_train = 20;
num_guard = 5;
P_fa = 0.1;
OS_detector = phased.CFARDetector('NumTrainingCells', 2 * num_train, ...
    'NumGuardCells', 2 * num_guard, 'ProbabilityFalseAlarm', P_fa);
OS_detector.Method = "OS";
OS_detector.Rank = num_train;
OS_detector.ThresholdOutputPort = true;

SOCA_detector = phased.CFARDetector('NumTrainingCells', 2 * num_train, ...
    'NumGuardCells', 2 * num_guard, 'ProbabilityFalseAlarm', P_fa);
SOCA_detector.Method = "SOCA";
SOCA_detector.ThresholdOutputPort = true;

% tracking parameters
P_D = 0.9;
P_G = 1e-6;
L_ite = 5;
Sigma_vec = 1;
T_ver = 1;
Max_hyp = 6;

% initialization
Tar_set = [];
Clu_set = [];
time_cnt = 0;
Tar_his = [];
cbf_n = {};
No_tar = 1;
Pro_ini = 1;

Tra_hyp = TraHypothesis(Tar_set, Clu_set, No_tar, Pro_ini, time_cnt, ...
    Tar_his, Sigma_vec, P_G, L_ite, T_ver, Max_hyp, 20);
Hypset_last = Tra_hyp;

tic;
for n_snap = 1 : Nsnap
    vec_n = Tra_mat(n_snap, :)';
    [~, mea_det] = findpeaks(vec_n);
    % det_diff = find(mea_det < 0);

    spec_rep = repmat(vec_n .^ 2, [3, 1]);
    [res_n, thres_n] = OS_detector(spec_rep, NB + 1 : 2 * NB);
    mea_now = intersect(mea_det, find(res_n));

    cbf_n{end + 1} = mea_now;

    % data association and update
    Hypset_snap = [];
    for h_idx = 1 : length(Hypset_last)
        [~, All_hyp] = Hypset_last(h_idx).SPA_track(mea_now);
        Hypset_snap = [Hypset_snap; All_hyp]; %#ok<AGROW> 
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

    if n_snap == 20
        1;
    end
    n_snap
end

Tra_hyp = Hypset_last(1);
for tar_idx = 1 : length(Tra_hyp.Tar_set)
    if(sum(Tra_hyp.Tar_set(tar_idx).his_Val) > Tra_hyp.Num_mote)
        Tra_hyp.Tra_his = [Tra_hyp.Tra_his; Tra_hyp.Tar_set(tar_idx)];
    end
end
toc;

figure; pcolor(Tra_mat);
colormap('jet'); shading flat; hold on; % clim([0, 3000]);
for s_idx = 1 : n_snap
    plot(cbf_n{s_idx}, ones(1, length(cbf_n{s_idx})) * s_idx, 'r.', ...
        'MarkerSize', 3);
end

state_his = {};
t_his = {};

LineColors = hsv(10);
labelOffset = 1;
figure; pcolor(Tra_mat);
colormap('gray'); shading flat;
hold on;
for tar_idx = 1 : length(Tra_hyp.Tra_his)
    plot(Tra_hyp.Tra_his(tar_idx).his_state, Tra_hyp.Tra_his(tar_idx).his_t, ...
        'LineWidth', 1, 'Color', LineColors(mod(tar_idx, 10) + 1, :));
    xLabelPos = mean(Tra_hyp.Tra_his(tar_idx).his_state) + labelOffset;
    yLabelPos = mean(Tra_hyp.Tra_his(tar_idx).his_t) + labelOffset;
    text(xLabelPos, yLabelPos, num2str(tar_idx), ...
        'Color', LineColors(mod(tar_idx, 10) + 1, :), ...
        'FontSize', 8);
    state_his{tar_idx} = Tra_hyp.Tra_his(tar_idx).his_state;
    t_his{tar_idx} = Tra_hyp.Tra_his(tar_idx).his_t;
end

% save('XuMenghuai_9_FZBTR6.mat', 'state_his', 't_his', 'cbf_n', 'Tra_mat');
% save('XuMenghuai_1_20250228_2x6B_1Tra', 'state_his', 't_his', 'cbf_n', 'Tra_mat');
% save('XuMenghuai_2_20250228_11_1.mat', 'state_his', 't_his', 'cbf_n', 'Tra_mat');




