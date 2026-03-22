clc; clear; close all; fclose all;
addpath(genpath("E:\\xmh_2025\\2505Tracking\\Matlab_tracking"));
lw = 1;
fsz = 12;
msz = 5;
clim_up = 100;
t_vec = linspace(0, 2 * pi, 360);

Nbeams = 453;
num_batch = 200;
len_snap = 512;

file_path = 'E:\\xmh_2025\\250922T\\data_bin\\';
fid_data = fopen([file_path, 'pdata.bin'], "r");

% detection initialization
N_tra = [500, 0];
N_gua = [50, 0];
alpha_dB = 11;
CA_det = CACFARdet(N_tra, N_gua);
Rist_min = 1e3;
P_G = 1e-6;
L_ite = 5;

Cov_mat = diag([100, 2]);

count_t = 0;
K_tar = 0;
CharaTar_set = [];
data_all = {};
while 1
    count_t = count_t + 1
    data_raw = fread(fid_data, [len_snap * num_batch * 2, Nbeams], 'float');
    if (size(data_raw, 2) < Nbeams), break; end
    data_mat = data_raw(1 : 2 : end, :) + 1j * data_raw(2 : 2 : end, :);
    data_mat(data_mat(:, 1) == 0, :) = [];
    if isempty(data_mat)
        continue;
    end
    mw_hfm = abs(data_mat);
    data_all{end + 1} = mw_hfm;
    h_fig1 = figure(1); clf(h_fig1); pcolor(mw_hfm);
    shading flat; colormap jet; clim([0 100]);
    N_range = size(data_mat, 1);
    Rist_max = N_range - N_tra(1);

    % detection
    mea_det = CA_det.detection2D(mw_hfm(Rist_min : Rist_max, :), alpha_dB);
    mea_det(:, 1) = mea_det(:, 1) + Rist_min - 1;

    % clustering
    [mea_clu, size_clu, res_clu] = cluster_link(mea_det, [1 30]);
    M_mea = size(mea_clu, 1);

    % plot the result of detection
    h_fig2 = figure(2); clf(h_fig2);
    h_pclr2 = pcolor(mw_hfm);
    shading flat; colormap(jet); hold on; clim([0 clim_up]);
    plot(mea_det(:, 2), mea_det(:, 1), 'w*', 'Markersize', msz);

    % plot the result of clustering
    for m_idx = 1 : M_mea
        mea_clun = mea_clu(m_idx, :);
        ellipse_pos = mea_clun' + [(size_clu(m_idx, 1)/2 + 20) * sin(t_vec);...
            (size_clu(m_idx, 2) / 2 + 0.5) * cos(t_vec)];
        figure(h_fig2); plot(ellipse_pos(2, :), ellipse_pos(1, :), ...
            'r-', 'LineWidth', lw);
    end
    figure(h_fig2); hold off;

    Measet_t = [];
    for m_idx = 1 : M_mea
        pointset_n = mea_det(res_clu == m_idx, :);
        amp_ind = sub2ind([N_range, Nbeams], pointset_n(:, 1), ...
            pointset_n(:, 2));
        ampset_n = mw_hfm(amp_ind);
        CharaMea_n = CharaMea(pointset_n, ampset_n);
        Measet_t = [Measet_t; CharaMea_n];
    end

    if (M_mea == 0) && (K_tar == 0)
        continue;
    elseif (M_mea > 0) && (K_tar == 0)
        for m_idx = 1 : M_mea
            CharaTar_k = CharaTarget(Measet_t(m_idx));
            CharaTar_set = [CharaTar_set; CharaTar_k];
        end
        K_tar = length(CharaTar_set);
    elseif (M_mea == 0) && (K_tar > 0)
    elseif (M_mea > 0) && (K_tar > 0)
        for k_idx = 1 : K_tar
            CharaTar_set(k_idx) = CharaTar_set(k_idx).predict();
        end
        beta_mat = charaset_associate(CharaTar_set, Measet_t, Cov_mat);
        beta_mat(beta_mat < P_G) = 0;
        udet_vec = P_G * ones(K_tar, 1);
        beta_mat = [beta_mat, udet_vec];

        xi_mat = [ones(K_tar, M_mea); 0.01 * ones(1, M_mea)];

        % SPA algorithm iterate
        phimat_itel = beta_mat(:, 1 : M_mea) ./ (beta_mat(:, 1 + M_mea) + eps);
        vmat_itel = zeros(K_tar, M_mea);
        for l_iter = 1 : L_ite
            vmat_itel = xi_mat(1 : K_tar, :) ./ (xi_mat(K_tar + 1, :) + ...
                sum(phimat_itel .* xi_mat(1 : K_tar, :), 1) - ...
                phimat_itel .* xi_mat(1 : K_tar, :));
            phimat_itel = beta_mat(:, 1 : M_mea) ./ (beta_mat(:, M_mea + 1) +...
                sum(beta_mat(:, 1 : M_mea) .* vmat_itel, 2) - beta_mat(:, ...
                1 : M_mea) .* vmat_itel);
        end
        % obtain the association probability
        Promat_a(:, 1 : M_mea) = beta_mat(:, 1 : M_mea) .* vmat_itel ./ ...
            (beta_mat(:, M_mea + 1) + sum(beta_mat(:, 1 : M_mea)...
            .* vmat_itel, 2));
        Promat_a(:, M_mea + 1) = beta_mat(:, M_mea + 1) ./ ...
            (beta_mat(:, M_mea + 1) + sum(beta_mat(:, 1 : M_mea)...
            .* vmat_itel, 2));
        [~, Pidx_max] = max(Promat_a, [], 2);

        K_ass = setdiff(1 : K_tar, find(Pidx_max == M_mea + 1))';
        M_ass = Pidx_max(Pidx_max <= M_mea);
        km_ass = [K_ass, M_ass];
        for km_idx = 1 : length(K_ass)
            k_ass = K_ass(km_idx);
            m_ass = M_ass(km_idx);
            if ismember(k_ass, K_ass(1 : km_idx - 1)), continue; end
            CharaTar_set(k_ass) = CharaTar_set(k_ass).update(Measet_t(m_ass));
        end
    end

    1;
end


