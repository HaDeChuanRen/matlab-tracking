clc; clear; close all; rng(1);

num_r = 20;
num_mc = 3000;
r_vec = 10 .^ linspace(-log10(4), log10(4), num_r)';
var_x = 1;
var_y = 1;
Cov_mat = diag([var_x, var_y]);

lambda_tar = 7;
kappa_cov = 10;
observe_cov = [var_x var_y] / kappa_cov;

P_Eu = zeros(num_r, num_mc);
P_Cov = zeros(num_r, num_mc);
P_Point2Point = zeros(num_r, num_mc);
P_Cha2Point = zeros(num_r, num_mc);
P_Cha2Cha = zeros(num_r, num_mc);

tic;
for mc_idx = 1 : num_mc
    h_waitbar = waitbar(mc_idx / num_mc);
    for r_idx = 1 : num_r
        r = r_vec(r_idx);
        d_cen = 5 * rand();
        theta_cen = 2 * pi * rand();
        point_cen = [d_cen * cos(theta_cen), d_cen * sin(theta_cen)];
        theta_dir = 2 * pi * rand();
        point_ref = point_cen + [r * sin(theta_dir), r * cos(theta_dir)];
        num_Points = poissrnd(lambda_tar) + 1;
        dis_points = [observe_cov .* randn(num_Points, 2); ...
            observe_cov .* randn(num_Points, 2) + [var_x * 1.5, 0]; ...
            observe_cov .* randn(num_Points, 2) - [var_x * 1.5, 0]];
        Set_Cen = point_cen + dis_points;
        Set_Ref = point_ref + dis_points;
        centroid_cen = mean(Set_Cen, 1);
        centroid_ref = mean(Set_Ref, 1);
        delta_cf = centroid_cen - centroid_ref;

        P_Eu(r_idx, mc_idx) = exp(-norm(delta_cf) ^ 2 / 2);
        P_Cov(r_idx, mc_idx) = exp(-(delta_cf / Cov_mat * delta_cf') / 2);
        Poicen_ChaMea = CharaMea(centroid_cen, 1);
        PoiRef_ChaMea = CharaMea(centroid_ref, 1);
        Poicen_ChaTar = CharaTarget(Poicen_ChaMea, 101, 1);
        % PoiRef_ChaTar = CharaTarget(PoiRef_ChaMea, 102, 1);
        P_Point2Point(r_idx, mc_idx) = chara_associate(Poicen_ChaTar, ...
            PoiRef_ChaMea, Cov_mat, [100; 100], 5);

        Setcen_ChaMea = CharaMea(Set_Cen, ones(num_Points * 3, 1));
        Setcen_ChaTar = CharaTarget(Setcen_ChaMea, 201, 1);
        Setref_ChaMea = CharaMea(Set_Ref, ones(num_Points * 3, 1));
        % Setref_ChaTar = CharaTarget(Setref_ChaMea, 202, 1);
        P_Cha2Point(r_idx, mc_idx) = chara_associate(Setcen_ChaTar, ...
            PoiRef_ChaMea,Cov_mat, [100; 100], 5);
        P_Cha2Cha(r_idx, mc_idx) = chara_associate(Setcen_ChaTar, ...
            Setref_ChaMea, Cov_mat / 5, [100; 100], 5);
    end
end
run_time = toc;
delete(h_waitbar);

if (run_time > 3600)
    % 获取当前时刻作为mat文件名的一部分
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    save(['save_mat/PvsR_' timestamp '.mat'], 'P_Eu', 'P_Cov', ...
        'P_Point2Point', 'P_Cha2Point','P_Cha2Cha', 'r_vec');
end

Pave_Eu = mean(P_Eu, 2);
Pave_Cov = mean(P_Cov, 2);
Pave_Point2Point = mean(P_Point2Point, 2);
Pave_Cha2Point = mean(P_Cha2Point, 2);
Pave_Cha2Cha = mean(P_Cha2Cha, 2);

lw = 1; msz = 10; fsz = 10;
figure; hold on; grid on;
semilogx(r_vec, Pave_Cov, 'bx--', 'LineWidth', lw);
semilogx(r_vec, Pave_Point2Point, 'r-', 'LineWidth', lw);
semilogx(r_vec, Pave_Cha2Point, 'm-', 'LineWidth', lw);
semilogx(r_vec, Pave_Cha2Cha, 'g-', 'LineWidth', lw);
legend('马氏距离', 'KL散度(点对点)', 'KL散度(点对扩展)', 'KL散度(扩展对扩展)');


