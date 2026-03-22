function beta_mat = charaset_associate(CharaTra_set, CharaMea_set, Cov_mat)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    K_tar = length(CharaTra_set);
    M_mea = length(CharaMea_set);
    beta_mat = zeros(K_tar, M_mea);

    for k_idx = 1 : K_tar
        for m_idx = 1 : M_mea
            beta_km = chara_associate(CharaTra_set(k_idx), ...
                CharaMea_set(m_idx), Cov_mat);
            beta_mat(k_idx, m_idx) = beta_km;
        end
    end

end


