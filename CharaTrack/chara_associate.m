function beta_ass = chara_associate(charaTar, charaMea, ...
    Cov_mat, N_grid, Sigma_rate)
% chara_associate 计算特征目标与特征观测的关联概率

    arguments
        charaTar (1,1) CharaTrack       % 输入的特征目标
        charaMea (1,1) CharaMea         % 输入的特征观测
        Cov_mat (:,:) double = []       % 关联协方差矩阵 Dim * Dim
        N_grid (:,1) double = []        % KL散度网格点数
        Sigma_rate (1,1) double = 3     % 关联方差扩大倍数
    end

    D_dim = charaTar.Dim;
    if isempty(Cov_mat), Cov_mat = eye(D_dim); end

    if isempty(N_grid), N_grid = 10 * ones(D_dim, 1);
    elseif (length(N_grid) < D_dim)
        N_grid = [N_grid; N_grid(end) * ones(D_dim - length(N_grid), 1)];
    end

    point_tar = charaTar.PointSet;
    amp_tar = charaTar.AmpSet / sum(charaTar.AmpSet);

    point_mea = charaMea.PointSet;
    amp_mea = charaMea.AmpSet / sum(charaMea.AmpSet);

    bound_tar = zeros(D_dim, 2);
    for d_idx = 1 : D_dim
        bound_tar(d_idx, 1) = floor(min(point_tar(:, d_idx))) - ...
            ceil(sqrt(Cov_mat(d_idx, d_idx))) * Sigma_rate;
        bound_tar(d_idx, 2) = ceil(max(point_tar(:, d_idx))) + ...
            ceil(sqrt(Cov_mat(d_idx, d_idx))) * Sigma_rate;
    end
    range_tar = bound_tar(:, 2) - bound_tar(:, 1);
    num_cells = prod(N_grid);
    grid_area = prod(range_tar ./ N_grid);

    ind_vec = zeros(D_dim, 1);
    KL_sum = 0;
    debug_prob = zeros(N_grid(1), N_grid(2));
    debug_ind = zeros(D_dim, 1);
    for c_idx = 1 : num_cells
        c_now = c_idx;
        for d_idx = 1 : D_dim
            debug_ind(d_idx) = mod(c_now, N_grid(d_idx)) + 1;
            ind_vec(d_idx) = ((mod(c_now, N_grid(d_idx)) + 0.5) / ...
                N_grid(d_idx)) * range_tar(d_idx) + bound_tar(d_idx, 1);
            c_now = floor(c_now / N_grid(d_idx));
        end
        probc_tar = amp_tar' * exp(-diag((point_tar - ind_vec') / Cov_mat * ...
            (point_tar' - ind_vec)) / 2) / ((2 * pi) ^ (D_dim / 2) * ...
            sqrt(det(Cov_mat)));
        probc_mea = amp_mea' * exp(-diag((point_mea - ind_vec') / Cov_mat * ...
            (point_mea' - ind_vec)) / 2) / ((2 * pi) ^ (D_dim / 2) * ...
            sqrt(det(Cov_mat))) + realmin;
        KL_sum = KL_sum + probc_tar * log(probc_tar / probc_mea) * grid_area;
        debug_prob(debug_ind(1), debug_ind(2)) = probc_tar;
    end
    ind_x = range_tar(1) * (0 : N_grid(1) - 1) / N_grid(1) + bound_tar(1, 1);
    ind_y = range_tar(2) * (0 : N_grid(2) - 1) / N_grid(2) + bound_tar(2, 1);
    [ind_X, ind_Y] = meshgrid(ind_x, ind_y);
    % figure; pcolor(ind_X, ind_Y, debug_prob'); shading flat; colormap jet; axis equal;
    % figure; plot(point_tar(:, 1), point_tar(:, 2), '*'); axis equal;
    beta_ass = exp(-KL_sum);
    1;
end