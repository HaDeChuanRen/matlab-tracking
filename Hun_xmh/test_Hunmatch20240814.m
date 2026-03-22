clc; clear; close all;
MC_times = 1000;
rng(1);
% cost_mat = [4 1 3 ; 2 0 5 ; 3 2 2 ];
cost_mat1 = eye(3);
cost_mat2 = [0 1 1 1; 1 0 1 1; 1 1 0 1; 1 1 1 0];

for mc_idx = 1 : MC_times
    n_mat = randi(8) + 2;
    m_mat = randi(8) + 2;
    cost_mat = randn(n_mat, m_mat);
    [assign_Hun, cost_Hun] = Hunmatch_xmh(cost_mat);
    [assign_Bao, cost_Bao] = Baomatch_xmh(cost_mat);
    if cost_Hun ~= cost_Bao
        dt_now = datetime("now");
        dt_now.Format = 'yyyyMMdd_HHmmss';
        file_name = ['costmat_', char(dt_now), '.mat'];
        save(file_name, "cost_mat");
        1;
    end
end