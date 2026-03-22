clc; clear; close all; rng(1);
lw = 1;
fsz = 12;
msz = 5;
t_vec = linspace(0, 2 * pi, 360);
LineColors = hsv(10);

%% targets location, velocity, appear time disapear time initialization
sample_time = 60;   % the total samples of tracking
T = 0.1;    % sample interval
sigma_w = [1e-6; 1e-6];
sigma_v = [0.001; 0.001];    % measurement noise variance
K_targets = 1;      % number of targets
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
% state transition matrix
lambda_tar = 10;
lambda_clu = 10;
observe_cov = [0.1 0.1];
No_tol = 1;

Amat = [1, T, 0, 0; 0, 1, 0, 0; 0, 0, 1, T; 0, 0, 0, 1];
Gamat_w = [T^2 / 2, 0; T, 0; 0, T^2; 0, T];

Xtrue_kt = Xori_mat + Gamat_w * (sqrt(sigma_w) .* randn(2, 1));
ground_all = zeros(4, sample_time, K_targets);

% Scene generating
dets_all = {};
for t_time = 1 : sample_time
    ground_all(1, t_time, :) = Xtrue_kt(1, :);
    ground_all(2, t_time, :) = Xtrue_kt(2, :);
    ground_all(3, t_time, :) = Xtrue_kt(3, :);
    ground_all(4, t_time, :) = Xtrue_kt(4, :);

    dets_t = zeros(0, 2);
    for k_idx = 1 : K_targets
        pointnum_kt = poissrnd(lambda_tar);
        tar_cen = [Xtrue_kt(1, k_idx), Xtrue_kt(3, k_idx)];
        det_kt = tar_cen + observe_cov .* randn(pointnum_kt, 2);
        dets_t = [dets_t; det_kt];
    end
    clunum_t = poissrnd(lambda_clu);
    xclu_t = xlim_l + (xlim_r - xlim_l) * rand(clunum_t, 1);
    yclu_t = ylim_d + (ylim_u - ylim_d) * rand(clunum_t, 1);
    clutters_t = [xclu_t, yclu_t];
    dets_t = [dets_t; clutters_t];
    dets_all{end + 1} = dets_t;
    % figure(1); clf; plot(dets_t(:, 1), dets_t(:, 2), '*');
    % grid on; xlim([xlim_l xlim_r]); ylim([ylim_d ylim_u]);
    % pause(0.3);

    Xtrue_kt = Amat * Xtrue_kt + Gamat_w * (sqrt(sigma_w) .* randn(2, 1));
end

%% Tracking
Cov_mat = diag([1, 1]);
P_G = 1e-6;
L_ite = 5;
CharaTar_set = [];
for t_time = 1 : sample_time
    dets_t = dets_all{t_time};
    det_dist = pdist(dets_t);
    det_link = linkage(det_dist, 'average');
    cluind_t = cluster(det_link, 'cutoff', 0.5, 'criterion', 'distance');

    ChaMeas_t = [];
    h_fig2 = figure(2); clf; plot(dets_t(:, 1), dets_t(:, 2), '*');
    grid on; hold on; xlim([xlim_l xlim_r]); ylim([ylim_d ylim_u]);
    title(num2str(t_time));
    for c_idx = 1 : max(cluind_t)
        dets_c = dets_t(cluind_t == c_idx, :);
        num_points = size(dets_c, 1);
        ChaMea_mt = CharaMea(dets_c, ones(num_points, 1));
        ChaMeas_t = [ChaMeas_t; ChaMea_mt];
        if num_points == 1
            mea_clun = dets_c;
        else
            mea_clun = mean(dets_c);
        end
        size_x = max(dets_c(:, 1)) - min(dets_c(:, 1)) + 0.4;
        size_y = max(dets_c(:, 2)) - min(dets_c(:, 2)) + 0.4;
        ellipse_pos = mea_clun' + [size_x / 2 * sin(t_vec); size_y / 2 * ...
            cos(t_vec)];
        figure(h_fig2); plot(ellipse_pos(1, :), ellipse_pos(2, :), ...
            'r-', 'LineWidth', lw);
    end

    M_mea = size(ChaMeas_t, 1);
    K_tar = size(CharaTar_set, 1);
    Promat_a = zeros(K_tar, M_mea + 1);
    Promat_b = zeros(K_tar + 1, M_mea);
    if (M_mea == 0) && (K_tar == 0)
        continue;
    elseif (M_mea > 0) && (K_tar == 0)
        for m_idx = 1 : M_mea
            CharaTar_k = CharaTarget(ChaMeas_t(m_idx), No_tol);
            No_tol = No_tol + 1;
            CharaTar_set = [CharaTar_set; CharaTar_k];
        end
        K_tar = length(CharaTar_set);
    elseif (M_mea == 0) && (K_tar > 0)
    elseif (M_mea > 0) && (K_tar > 0)
        for k_idx = 1 : K_tar
            CharaTar_set(k_idx) = CharaTar_set(k_idx).predict();
        end
        beta_mat = charaset_associate(CharaTar_set, ChaMeas_t, Cov_mat);
        beta_mat(beta_mat < P_G) = 0;
        udet_vec = P_G * ones(K_tar, 1);
        beta_mat = [beta_mat, udet_vec];

        % SPA algorithm iterate
        % xi_mat = [ones(K_tar, M_mea); 0.01 * ones(1, M_mea)];
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
        % Promat_a(:, 1 : M_mea) = beta_mat(:, 1 : M_mea) .* vmat_itel ./ ...
        %     (beta_mat(:, M_mea + 1) + sum(beta_mat(:, 1 : M_mea)...
        %     .* vmat_itel, 2));
        % Promat_a(:, M_mea + 1) = beta_mat(:, M_mea + 1) ./ ...
        %     (beta_mat(:, M_mea + 1) + sum(beta_mat(:, 1 : M_mea)...
        %     .* vmat_itel, 2));
        % [~, Pidx_max] = max(Promat_a, [], 2);

        [~, Pidx_max] = max(beta_mat, [], 2);
        cost_beta = 1 - beta_mat(:, 1 : end - 1);
        assign_mat = Hunmatch_xmh(cost_beta);
        for a_idx = 1 : size(assign_mat, 1)
            if (Pidx_max(assign_mat(a_idx, 1)) == M_mea + 1)
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
            TarObj_n = CharaTarget(ChaMeas_t(m_idx), No_tol);
            No_tol = No_tol + 1;
            CharaTar_set = [CharaTar_set; TarObj_n];
        end

        % banish the undetected targets
        K_ban = [];
        for k_idx = 1 : K_tar
            % if ismember(k_idx, K_ass), continue; end
            if CharaTar_set(k_idx).P_s <= 0
                K_ban(end + 1) = k_idx;
            end
        end
        CharaTar_set(K_ban) = [];
        K_tar = length(CharaTar_set);
        K_show = [];
        for k_idx = 1 : K_tar
            if CharaTar_set(k_idx).P_s >= 0.5
                K_show(end + 1) = k_idx;
            end
        end

        for ks_idx = K_show
            figure(h_fig2); plot(CharaTar_set(ks_idx).his_peak(:, 1),...
                CharaTar_set(ks_idx).his_peak(:, 2), 'LineWidth', 1, 'Color',...
                LineColors(mod(CharaTar_set(ks_idx).ID_tar, 10) + 1, :));
        end
    end
    % pause(1);
    if t_time == 39
        1;
    end
end





