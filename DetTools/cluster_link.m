function [mea_clu, size_clu, res_clu] = cluster_link(mea_det, expand_vec)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[N_det, D_dim] = size(mea_det);
if N_det == 0
    mea_clu = [];
    size_clu = [];
    res_clu = [];
    return;
elseif N_det == 1
    mea_clu = mea_det;
    size_clu = zeros(1, D_dim);
    res_clu = 1;
    return;
end

mea_temp = mea_det .* expand_vec;
Dis_tar = pdist(mea_temp);

Fre_link = linkage(Dis_tar, 'single'); %figure; dendrogram(Fre_link)
res_clu = cluster(Fre_link, 'Cutoff', 100, 'Criterion', 'distance');
N_max = max(res_clu);

mea_clu = zeros(N_max, D_dim);
size_clu = zeros(N_max, D_dim);

for n_idx = 1 : N_max
    idx_clun = find(res_clu == n_idx);
    if isempty(idx_clun)
        continue;
    end
    for d_idx = 1 : D_dim
        size_d = max(mea_det(idx_clun, d_idx)) - min(mea_det(idx_clun, d_idx));
        size_clu(n_idx, d_idx) = size_d;
    end

    mea_clun = mean(mea_det(idx_clun, :), 1);
    mea_clu(n_idx, :) = mea_clun;
end


end