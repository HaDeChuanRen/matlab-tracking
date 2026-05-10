% 构建跟踪算法测试场景并可视化展示
% used to build the test scene for tracking algorithm and show it visibly
clear; clc; close all;
rng(1);
addpath(genpath('C:\study\Codes\Matlab\Tracking'));

% plot parameters
lw = 1;
fsz = 12;
msz = 5;
t_vec = linspace(0, 2 * pi, 360);
LineColors = hsv(10);

%% targets location, velocity, appear time disapear time initialization
sample_time = 60;   % the total samples of tracking
T = 0.1;    % sample interval
K_targets = 6;      % number of targets
sigma_w = [1e-6; 1e-6];     % noise variance of localization
sigma_v = [0.001; 0.001];   % noise variance of velocity

xrange = 10;
xmin = -5;
xPosition_true0 = xrange * rand(K_targets, 1) + xmin;
vxrange = 2;
vxmin = -1;
xVelocity_true0 = vxrange * rand(K_targets, 1) + vxmin;
yrange = 16;
ymin = -8;
yPosition_true0 = yrange * rand(K_targets, 1) + ymin;
vyrange = 4;
vymin = -2;
yVelocity_true0 = vyrange * rand(K_targets, 1) + vymin;
Xori_mat = [xPosition_true0'; xVelocity_true0'; yPosition_true0'; ...
    yVelocity_true0'];

xlim_l = xmin + sample_time * min([vxmin, 0]) * T;
xlim_r = xmin + xrange + sample_time * max([vxrange + vxmin, 0]) * T;
ylim_d = ymin + sample_time * min([vymin, 0]) * T;
ylim_u = ymin + yrange + sample_time * max([vyrange + vymin, 0]) * T;

%% create targets and clutters in the scene
lambda_tar = 4;
lambda_clu = 10;
observe_cov = [0.1 0.1];
No_tol = 1;

% state transition matrix
Amat = [1, T, 0, 0; 0, 1, 0, 0; 0, 0, 1, T; 0, 0, 0, 1];
Gamat_w = [T^2 / 2, 0; T, 0; 0, T^2; 0, T];

% initialize the baseTrack for ground truth of targets
Xtrue_kt = Xori_mat + Gamat_w * (sqrt(sigma_w) .* randn(2, 1));
for k_idx = 1 : K_targets
    allTrueTrack(k_idx) = BaseClass.BaseTrack(k_idx, 2, 1, ...
        [Xtrue_kt(1, k_idx), Xtrue_kt(3, k_idx)]);
end

% Scene generating
dets_all = {};
allShapes = ['|', '-', '+', 'd', 's', '<', '>', '^', 'v', 'h'];
numShapes = length(allShapes);
TarShapes = allShapes(ceil(numShapes * rand(K_targets, 1)));
for t_time = 1 : sample_time
    % use basetrack to store the ground truth of targets
    if t_time > 1
        Xtrue_kt = Amat * Xtrue_kt + Gamat_w * (sqrt(sigma_w) .* randn(2, 1));
        for k_idx = 1 : K_targets
            allTrueTrack(k_idx) = allTrueTrack(k_idx).base_add(...
                [Xtrue_kt(1, k_idx), Xtrue_kt(3, k_idx)], t_time);
        end
    end

    dets_t = zeros(0, 2);
    for k_idx = 1 : K_targets
        tar_cen = [Xtrue_kt(1, k_idx), Xtrue_kt(3, k_idx)];
        % pointnum_kt = poissrnd(lambda_tar);
        % det_kt = tar_cen + observe_cov .* randn(pointnum_kt, 2);
        det_kt = charaDets2D(TarShapes(k_idx), tar_cen, observe_cov);
        dets_t = [dets_t; det_kt];
    end
    clunum_t = poissrnd(lambda_clu);
    xclu_t = xlim_l + (xlim_r - xlim_l) * rand(clunum_t, 1);
    yclu_t = ylim_d + (ylim_u - ylim_d) * rand(clunum_t, 1);
    clutters_t = [xclu_t, yclu_t];
    dets_t = [dets_t; clutters_t];
    dets_all{end + 1} = dets_t;
end

%% Tracking and Plotting
% initializing
Cov_mat = diag([1, 1]);
P_G = 1e-6;
L_ite = 5;
CharaTar_set = [];
allShowNo = [];

tic;
h_fig1 = figure(1);
xlabel('xaxis(m)'); ylabel('yaxis(m)');
numColors = 12;

for t_time = 1 : sample_time
    waitbar(t_time / sample_time);

    % obtian measurements of this moment
    dets_t = dets_all{t_time};

    % plot real tracks and measurements
    clf(h_fig1);
    h_fig1; hold on;
    plot(dets_t(:, 1), dets_t(:, 2), 'r*', 'MarkerSize', msz);
    for k_idx = 1 : K_targets
        plot(allTrueTrack(k_idx).TrackInfo(1:t_time, 1), ...
            allTrueTrack(k_idx).TrackInfo(1:t_time, 2), 'k--','LineWidth', lw);
    end
    h_fig1; hold off; grid on;
    xlim([xlim_l xlim_r]); ylim([ylim_d ylim_u]);

    % clustering
    det_dist = pdist(dets_t);
    det_link = linkage(det_dist, 'average');
    cluind_t = cluster(det_link, 'cutoff', 0.5, 'criterion', 'distance');
    ChaMeas_t = [];
    allMeaColors = hsv(max(cluind_t));
    for c_idx = 1 : max(cluind_t)
        dets_c = dets_t(cluind_t == c_idx, :);
        num_points = size(dets_c, 1);
        ChaMea_mt = CharaMea(dets_c, ones(num_points, 1));
        ChaMeas_t = [ChaMeas_t; ChaMea_mt];
        h_fig1; hold on;
        h_fig1; plot(dets_c(:, 1), dets_c(:, 2), '*', 'MarkerSize', msz, ...
            'Color', allMeaColors(c_idx, :));
        h_fig1; hold off;
    end

    % considering
    M_mea = size(ChaMeas_t, 1);
    K_tar = size(CharaTar_set, 1);
    if (M_mea == 0) && (K_tar == 0)
        continue;
    elseif (M_mea > 0) && (K_tar == 0)
        for m_idx = 1 : M_mea
            CharaTar_k = CharaTrack(ChaMeas_t(m_idx), No_tol);
            No_tol = No_tol + 1;
            CharaTar_set = [CharaTar_set; CharaTar_k];
        end
        K_tar = length(CharaTar_set);
    elseif (M_mea == 0) && (K_tar > 0)
        % prediction
        for k_idx = 1 : K_tar
            CharaTar_set(k_idx) = CharaTar_set(k_idx).predict();
        end
    elseif (M_mea > 0) && (K_tar > 0)
        % prediction
        for k_idx = 1 : K_tar
            CharaTar_set(k_idx) = CharaTar_set(k_idx).predict();
        end

        % obtain probabilities of association
        beta_mat = zeros(K_tar, M_mea);
        for k_idx = 1 : K_tar
            for m_idx = 1 : M_mea
                beta_km = chara_associate(CharaTar_set(k_idx), ...
                    ChaMeas_t(m_idx), Cov_mat, 20, 5);
                beta_mat(k_idx, m_idx) = beta_km;
            end
        end
        beta_mat(beta_mat < P_G) = 0;
        udet_vec = P_G * ones(K_tar, 1);
        beta_mat = [beta_mat, udet_vec];
        % figure(11); pcolor(beta_mat); colormap gray; shading flat;

        % SPA algorithm iterate
        % xi_mat = [ones(K_tar, M_mea); ones(1, M_mea)];
        % phimat_itel = beta_mat(:, 1 : M_mea) ./ (beta_mat(:, 1 + M_mea) + eps);
        % vmat_itel = zeros(K_tar, M_mea);
        % for l_iter = 1 : L_ite
        %     vmat_itel = xi_mat(1 : K_tar, :) ./ (xi_mat(K_tar + 1, :) + ...
        %         sum(phimat_itel .* xi_mat(1 : K_tar, :), 1) - ...
        %         phimat_itel .* xi_mat(1 : K_tar, :));
        %     phimat_itel = beta_mat(:, 1 : M_mea) ./ (beta_mat(:, M_mea + 1) +...
        %         sum(beta_mat(:, 1 : M_mea) .* vmat_itel, 2) - beta_mat(:, ...
        %         1 : M_mea) .* vmat_itel);
        % end
        % % obtain the association probability
        % Promat_a = zeros(K_tar, M_mea + 1);
        % Promat_b = zeros(K_tar + 1, M_mea);
        % Promat_a(:, 1 : M_mea) = beta_mat(:, 1 : M_mea) .* vmat_itel ./ ...
        %     (beta_mat(:, M_mea + 1) + sum(beta_mat(:, 1 : M_mea)...
        %     .* vmat_itel, 2));
        % Promat_a(:, M_mea + 1) = beta_mat(:, M_mea + 1) ./ ...
        %     (beta_mat(:, M_mea + 1) + sum(beta_mat(:, 1 : M_mea)...
        %     .* vmat_itel, 2));
        % % figure(12); pcolor(Promat_a); colormap gray; shading flat;
        % [~, Pidx_max] = max(Promat_a, [], 2);

        [~, Pidx_max] = max(beta_mat, [], 2);
        cost_beta = - log(beta_mat(:, 1 : end - 1) + eps);
        assign_mat = Hungarian(cost_beta);
        for a_idx = 1 : size(assign_mat, 1)
            if (Pidx_max(assign_mat(a_idx, 1)) == M_mea + 1 || ...
                    beta_mat(assign_mat(a_idx, 1), assign_mat(a_idx, 2)) < P_G)
                assign_mat(a_idx, 2) = M_mea + 1;
            end
        end
        assign_mat(assign_mat(:, 2) == M_mea + 1, :) = [];
        K_ass = assign_mat(:, 1);
        M_ass = assign_mat(:, 2);
        for km_idx = 1 : length(K_ass)
            k_ass = K_ass(km_idx);
            m_ass = M_ass(km_idx);
            if ismember(k_ass, K_ass(1 : km_idx - 1)), continue; end
            CharaTar_set(k_ass) = CharaTar_set(k_ass).update(ChaMeas_t(m_ass));
        end

        % create new targets by measurements
        for m_idx = 1 : M_mea
            if ismember(m_idx, M_ass), continue; end
            TarObj_n = CharaTrack(ChaMeas_t(m_idx), No_tol);
            No_tol = No_tol + 1;
            CharaTar_set = [CharaTar_set; TarObj_n];
        end

        % banish the undetected targets
        K_ban = [];
        for k_idx = 1 : K_tar
            % if ismember(k_idx, K_ass), continue; end
            if CharaTar_set(k_idx).Ps <= 0
                K_ban(end + 1) = k_idx;
            end
        end
        CharaTar_set(K_ban) = [];
        K_tar = length(CharaTar_set);
        K_show = [];
        for k_idx = 1 : K_tar
            if CharaTar_set(k_idx).Ps >= 0.5
                K_show(end + 1) = k_idx;
            end
        end

        for ks_idx = K_show
            [~, locB] = ismember(CharaTar_set(ks_idx).TrackID, allShowNo);
            if locB > 0
                allShowTrack(locB).base_add(...
                    CharaTar_set(ks_idx).TrackInfo(end, :), ...
                    CharaTar_set(ks_idx).MomentInfo(end));
            else
                allShowNo(end + 1) = CharaTar_set(ks_idx).TrackID;
                numShow = length(allShowNo);
                allShowTrack(numShow) = BaseClass.BaseTrack(...
                    CharaTar_set(ks_idx).TrackID, ...
                    CharaTar_set(ks_idx).Dim, ...
                    CharaTar_set(ks_idx).MomentInfo, ...
                    CharaTar_set(ks_idx).TrackInfo);
            end
        end
    end
    1;
end
toc;
