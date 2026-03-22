function [ospa_pd, ospa_pfa, match_res] = OSPA_cal(ground_all, traj_all, time_stamp, ...
    C_thr, a_pun)
% OSPA_cal: calculate the optimal subpattern assignment for the tracking result
% Input:
% Output:

% Last update: 2025/8/21
% 2025/8/21 update details
% 1. (Todo) write the function explanations and the comments
% 2. (Todo) use the label of targets rewrite the calculation of OSPA
% 3. (Todo) revise the function by new port

% if ~exist('Delta_thr', 'var'), Delta_thr = 10;
% elseif isempty(Delta_thr), Delta_thr = 10; end

if ~exist('C_thr', 'var'), C_thr = 10;
elseif isempty(C_thr), C_thr = 10; end

if ~exist('a_pun', 'var'), a_pun = 1;
elseif isempty(a_pun), a_pun = 1; end

% state dimension and number of targets
k_targets = length(ground_all);
d_dim = size(ground_all{1, 1}, 2);

% estimated number of targets
k_est = length(traj_all);

% length of time stamp
num_frame = length(time_stamp);


% calculate the distance between each estimated trajectory and the truth
% trajectory
gepur_ten = zeros(k_targets, k_est, num_frame);
for k_idx = 1 : k_targets
    for k_jdx = 1 : k_est
        G_mati = ground_all{k_idx, 1};
        GE_mati = nan(num_frame, d_dim);
        Gstart_idx = find(time_stamp == ground_all{k_idx, 2}(1));
        Gend_idx = find(time_stamp == ground_all{k_idx, 2}(end));
        GE_mati(Gstart_idx : Gend_idx, :) = G_mati;

        T_matj = traj_all{k_jdx};
        TE_matj = nan(num_frame, d_dim);
        Tstart_idx = find(time_stamp == traj_all{k_jdx, 2}(1));
        Tend_idx = find(time_stamp == traj_all{k_jdx, 2}(end));
        TE_matj(Tstart_idx : Tend_idx, :) = T_matj;

        GTdel_mat = GE_mati - TE_matj;
        norm_vec = vecnorm(GTdel_mat, 2, 2);
        norm_vec(norm_vec > C_thr) = C_thr;
        norm_vec(isnan(norm_vec)) = C_thr;
        gepur_ten(k_idx, k_jdx, :) = norm_vec;
    end
end

gepur_mat = sum(gepur_ten, 3);
[match_res, ~] = Hunmatch_xmh(gepur_mat);

% calculate OSPA
ospa_dist = zeros(num_frame, 1);
for k_idx = 1 : k_targets
    ospa_kt = zeros(num_frame, 1);
    Gstart_idx = find(time_stamp == ground_all{k_idx, 2}(1));
    Gend_idx = find(time_stamp == ground_all{k_idx, 2}(end));
    k_midx = find(match_res(:, 1) == k_idx, 1);
    if isempty(k_midx)
        ospa_kt(Gstart_idx : Gend_idx) = C_thr;
    end
    k_jdx = match_res(k_midx, 2);
    ospa_kt(Gstart_idx : Gend_idx) = gepur_ten(k_idx, k_jdx, ...
        Gstart_idx : Gend_idx);
    ospa_dist = ospa_dist + ospa_kt;
end

ospa_pd = ospa_dist / k_targets;

for f_idx = 1 : k_est
    if ismember(f_idx, match_res(:, 2)), continue; end
    ospa_ft = zeros(num_frame, 1);
    Tstart_idx = find(time_stamp == traj_all{f_idx, 2}(1));
    Tend_idx = find(time_stamp == traj_all{f_idx, 2}(end));
    ospa_ft(Tstart_idx : Tend_idx) = C_thr;
    ospa_dist = ospa_dist + ospa_ft;
end
ospa_pfa = ospa_dist / k_est;

end