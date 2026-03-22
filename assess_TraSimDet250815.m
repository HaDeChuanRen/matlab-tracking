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
glmb_plat = GLMBPlat();

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

    angle_det = theta_vec(mea_det(:, 1));
    mea_pos = [mea_det(:, 2) .* cos(angle_det), ...
        mea_det(:, 2) .* sin(angle_det)];
    % mea_pos = [tar_pos; clu_pos];
    mea_col{end + 1} = mea_pos; %#ok<SAGROW>

    % data association and update
    [SPAMHT_plat, traj_spa, ~] = SPAMHT_plat.track(mea_pos);
    % glmb_plat = glmb_plat.track(mea_pos, n_idx);

    1;
    % pause(0.5);
end
toc;

ground_all = cell(k_targets, 2);
for k_idx = 1 : k_targets
    ground_all{k_idx, 1} = squeeze(tar_track(k_idx, :, :))';
    ground_all{k_idx, 2} = (1 : num_frame)';
end

h_fig2 = figure(2); clf(h_fig2); hold on; grid on;
for k_idx = 1 : k_targets
    plot(ground_all{k_idx, 1}(:, 1), ground_all{k_idx, 1}(:, 2), ...
        'k--', 'LineWidth', lw);
end

for n_idx = 1 : num_frame
    plot(mea_col{n_idx}(:, 1), mea_col{n_idx}(:, 2), 'r*');
end

LineColors = hsv(10);
for t_idx = 1 : length(traj_spa)
    plot(traj_spa{t_idx}(:, 1), traj_spa{t_idx}(:, 2), 'LineWidth', 1, ...
        'Color', LineColors(mod(t_idx, 10) + 1, :));
end
h_fig2; hold off;

time_stamp = 1 : num_frame;
C_thr = sqrt(x_cluedge ^ 2 + y_cluedge ^ 2);
[ospaPd_spa, ospaPfa_spa, ~] = OSPA_cal(ground_all, traj_spa, time_stamp, ...
    C_thr);

h_fig3 = figure(3); clf(h_fig3); hold on; grid on;
plot(time_stamp, ospaPd_spa, 'LineWidth', lw);

h_fig4 = figure(4); clf(h_fig4); hold on; grid on;
plot(time_stamp, ospaPfa_spa, 'LineWidth', lw);
