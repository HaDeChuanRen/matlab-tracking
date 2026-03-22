function [Tar_set, Promat_a, km_ass, K_show, Tar_ban] = SPA_Track(Tar_set, Mea_now, ...
    time_cnt, P_G, insig1, insig2, L_ite)
% last update: 2025/2/21

% default values of parameters
if ~exist('P_G', 'var'), P_G = 1e-5;
elseif isempty(P_G), P_G = 1e-5; end
if ~exist('insig1', 'var'), insig1 = 100;
elseif isempty(insig1), insig1 = 100; end
if ~exist('insig2', 'var'), insig2 = 5;
elseif isempty(insig2), insig2 = 5; end
if ~exist('L_ite', 'var'), L_ite = 5;
elseif isempty(L_ite), L_ite = 5; end

% data association and update
% initialization
K_tar = length(Tar_set);
M_mea = size(Mea_now, 1);
Promat_a = zeros(K_tar, M_mea + 1);
Promat_b = zeros(K_tar + 1, M_mea);
km_ass = [];
K_show = [];
if K_tar == 0
    No_tol = 1;
    for m_idx = 1 : M_mea
        Tar_m = syhTarget(No_tol, Mea_now(m_idx, 1), Mea_now(m_idx, 2), ...
            time_cnt);
        Tar_set = [Tar_set; Tar_m]; %#ok<*AGROW> 
        No_tol = No_tol + 1;
    end
elseif K_tar > 0
    No_tol = Tar_set(end).No_tar + 1;
    tar_range = zeros(K_tar, 1);
    tar_angle = zeros(K_tar, 1);
    for k_idx = 1 : K_tar
        [Tar_set(k_idx), range_k, angle_k] = ...
            Tar_set(k_idx).state_transform(time_cnt);
        tar_range(k_idx) = range_k;
        tar_angle(k_idx) = angle_k;
    end
    if M_mea == 0
        for k_idx = 1 : K_tar
            if Tar_set(k_idx).P_s >= 0.6
                K_show(end + 1) = k_idx;
            end
        end
        return;
    end
    cost_1 = abs(tar_range - Mea_now(:, 1)');
    cost_2 = abs(tar_angle - Mea_now(:, 2)');
    Smat = cost_1 .^ 2 / (insig1 ^ 2) + ...
        cost_2 .^ 2 / (insig2 ^ 2);
    betamat = exp(-Smat);
    betamat(betamat < P_G) = 0;
    udet_vec = P_G * ones(K_tar, 1);
    betamat = [betamat, udet_vec];
    ximat = [ones(K_tar, M_mea); 0.01 * ones(1, M_mea)];

    % SPA algorithm iterate
    phimat_itel = betamat(:, 1 : M_mea) ./ (betamat(:, 1 + M_mea) + eps);
    vmat_itel = zeros(K_tar, M_mea);
    for l_iter = 1 : L_ite
        vmat_itel = ximat(1 : K_tar, :) ./ (ximat(K_tar + 1, :) + ...
            sum(phimat_itel .* ximat(1 : K_tar, :), 1) - ...
            phimat_itel .* ximat(1 : K_tar, :));
        phimat_itel = betamat(:, 1 : M_mea) ./ (betamat(:, M_mea + 1) + ...
            sum(betamat(:, 1 : M_mea) .* vmat_itel, 2) - betamat(:, ...
            1 : M_mea) .* vmat_itel);
    end

    % obtain the association probability
    Promat_a(:, 1 : M_mea) = betamat(:, 1 : M_mea) .* vmat_itel ./ ...
        (betamat(:, M_mea + 1) + sum(betamat(:, 1 : M_mea)...
        .* vmat_itel, 2));
    Promat_a(:, M_mea + 1) = betamat(:, M_mea + 1) ./ ...
        (betamat(:, M_mea + 1) + sum(betamat(:, 1 : M_mea)...
        .* vmat_itel, 2));

    Promat_b(1 : K_tar, :) = ximat(1 : K_tar, :) .* phimat_itel ./ ...
        (ximat(K_tar + 1, :) + sum(ximat(1 : K_tar, :)...
        .* phimat_itel, 1));
    Promat_b(K_tar + 1, :) = ximat(K_tar + 1, :) ./ ...
        (ximat(K_tar + 1, :) + sum(ximat(1 : K_tar, :)...
        .* phimat_itel, 1));

    [~, Pidx_max] = max(Promat_a, [], 2);
    K_ass = setdiff(1 : K_tar, find(Pidx_max == M_mea + 1))';
    M_ass = Pidx_max(Pidx_max <= M_mea);
    km_ass = [K_ass, M_ass];

    for km_idx = 1 : length(K_ass)
        k_ass = K_ass(km_idx);
        m_ass = M_ass(km_idx);
        if ismember(k_ass, K_ass(1 : km_idx - 1)), continue; end
        Tar_set(k_ass) = Tar_set(k_ass).Update(Mea_now(m_ass, 1), ...
            Mea_now(m_ass, 2));
    end

    % create new targets by measurements
    % Mind_new = setdiff(1 : N_clu, ass_mat(:, 2));
    for m_idx = 1 : M_mea
        if ismember(m_idx, M_ass), continue; end 
        TarObj_n = syhTarget(No_tol, Mea_now(m_idx, 1), Mea_now(m_idx, 2), ...
            time_cnt);
        No_tol = No_tol + 1;
        Tar_set = [Tar_set; TarObj_n];
    end

    % banish the undetected targets
    K_ban = [];
    for k_idx = 1 : K_tar
        % if ismember(k_idx, K_ass), continue; end
        if Tar_set(k_idx).P_s <= 0
            K_ban(end + 1) = k_idx;
        end
    end
%     Tar_ban = Tar_set(K_ban);
    Tar_set(K_ban) = [];
    K_tar = length(Tar_set);
    for k_idx = 1 : K_tar
        if Tar_set(k_idx).P_s >= 0.5
            K_show(end + 1) = k_idx;
        end
    end
end

end