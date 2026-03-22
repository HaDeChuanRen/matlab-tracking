classdef DimTarget
    %   last date: 2025/3/17
    %   
    
    properties
        No_tar      % 目标编号（ID） int
        t_cnt       % 当前时刻 double 或 datetime
        his_t       % 时刻历史 T_last * 1 double 或 datetime
        T_ver       % 追踪时间间隔 double
        N_last      % 目标持续周期数 int
        D_dim       % 目标维度 int

        State_vec   % 目标状态 D_dim * 1 double
        Vel_vec     % 目标状态变化率 D_dim * 1 double
        his_state   % 状态历史 T_last * Dim double
        Sigma_vec   % 目标测量方差 D_dim * 1 double
        Cov_mat     % 目标状态协方差 (2 * D_dim) * (2 * D_dim) double

        P_s         % 存活概率 double
        P_d         % 检测概率 double
        his_Ps      % 存活概率历史 T_last * 1 double
        his_Val     % 有效关联历史 T_last * 1 int
    end
    
    methods
        function obj = DimTarget(inNo, inMoment, inState, inSigma, inPs, ...
            inInterval, inPd)
            if ~exist('inPs', 'var'), inPs = 0.3;
            elseif isempty(inPs), inPs = 0.3; end
            if ~exist('inInterval', 'var'), inInterval = 1;
            elseif isempty(inInterval), inInterval = 1; end
            if ~exist('inPd', 'var'), inPd = 0.98;
            elseif isempty(inPd), inPd = 0.98; end

            obj.No_tar = inNo;
            obj.t_cnt = inMoment;
            obj.his_t = obj.t_cnt ;
            obj.T_ver = inInterval;
            obj.N_last = 1;

            obj.State_vec = inState;
            obj.D_dim = size(inState, 1);
            obj.Vel_vec = zeros(obj.D_dim, 1);
            obj.Sigma_vec = inSigma;
            obj.Cov_mat = diag([inSigma; inSigma]);
            obj.his_state = obj.State_vec';

            obj.P_s = inPs;
            obj.his_Ps = obj.P_s;
            obj.P_d = inPd;
            obj.his_Val = 1;
        end

        function [obj, State_now] = Predict(obj, inMoment, inSigma)

            if exist('inMoment', 'var')
                obj.t_cnt = inMoment;
            else
                obj.t_cnt = obj.t_cnt + obj.T_ver;
            end
            obj.his_t = [obj.his_t; obj.t_cnt];
            if exist('inSigma', 'var'), obj.Sigma_vec = inSigma; end

            Qmat = eye(2 * obj.D_dim);
            Tmat = [eye(obj.D_dim), obj.T_ver * eye(obj.D_dim); zeros(obj.D_dim), ...
                eye(obj.D_dim)];

            obj.State_vec = obj.State_vec + obj.Vel_vec * obj.T_ver;
            obj.his_state = [obj.his_state; obj.State_vec'];
            obj.N_last = obj.N_last + 1;
            obj.his_Val = [obj.his_Val; 0];
            obj.his_Ps = [obj.his_Ps; obj.P_s];
            obj.Cov_mat = Tmat * obj.Cov_mat * Tmat' + Qmat;
            State_now = obj.State_vec;
        end

        % function obj = softUpdate(obj, Zmat, Ass_vec, inPd)
        %     if ~exist('inPd', 'var'), inPd = 0.98;
        %     elseif isempty(inPd), inPd = 0.98; end
        %     obj.P_d = inPd;

        %     beta_vec = Ass_vec(1 : end - 1);
        %     beta_0 = Ass_vec(end);

        %     if (beta_0 < 0.5)
        %         obj.Val_vec(end) = 1;
        %     end

        %     obj.his_Ps = [obj.his_Ps; obj.P_s];
        %     obj.P_s = 1 - (1 - obj.P_s) * beta_0 / (1 - obj.P_d * obj.P_s);
        %     if (obj.P_s >= 1), obj.P_s = 1 - 1e-5; end

        %     hmat = [1; 0];
        %     Ervec_x = Zmat(:, 1) - obj.state_x(1);
        %     Ermean_x = Ervec_x' * beta_vec;
        %     Kmat_x = obj.cov_x * hmat / (hmat' * obj.cov_x * hmat + ...
        %         obj.sigma_x);
        %     post_x = obj.state_x + Kmat_x * Ermean_x;
        %     Ptilde_x = Kmat_x * (Ervec_x' * diag(beta_vec) * Ervec_x - ...
        %         Ermean_x^2) * Kmat_x';
        %     obj.cov_x = obj.cov_x - (1 - beta_0) * Kmat_x * (hmat' * ...
        %         obj.cov_x * hmat + obj.sigma_x) * Kmat_x' + Ptilde_x;
        %     obj.state_x(1) = Ass_vec' * [Zmat(:, 1); obj.state_x(1)];
        %     obj.state_x(2) = post_x(2);

        %     Ervec_y = Zmat(:, 2) - obj.state_y(1);
        %     Ermean_y = Ervec_y' * beta_vec;
        %     Kmat_y = obj.cov_y * hmat / (hmat' * obj.cov_y * hmat + ...
        %         obj.sigma_y);
        %     post_y = obj.state_y + Kmat_y * Ermean_y;
        %     Ptilde_y = Kmat_y * (Ervec_y' * diag(beta_vec) * Ervec_y - ...
        %         Ermean_y^2) * Kmat_y';
        %     obj.cov_y = obj.cov_y - (1 - beta_0) * Kmat_y * (hmat' * ...
        %         obj.cov_y * hmat + obj.sigma_y) * Kmat_y' + Ptilde_y;
        %     obj.state_y(1) = Ass_vec' * [Zmat(:, 2); obj.state_y(1)];
        %     obj.state_y(2) = post_y(2);
        % end

        function obj = hardUpdate(obj, Z_ass, Ass_pro, inPd)
            if ~exist('inPd', 'var'), inPd = 0.98;
            elseif isempty(inPd), inPd = 0.98; end
            obj.P_d = inPd;

            if (Ass_pro > 0.5)
                obj.his_Val(end) = 1;
            end

            obj.P_s = 1 - (1 - obj.P_s) * (1 - Ass_pro) / (1 - obj.P_d * ...
                obj.P_s);
            if (obj.P_s >= 1), obj.P_s = 1 - 1e-12; end
            obj.his_Ps(end) = obj.P_s;

            if ~isempty(Z_ass)
                Hmat = [eye(obj.D_dim), zeros(obj.D_dim)];
                Kmat = (obj.Cov_mat * Hmat') / (diag(obj.Sigma_vec) + Hmat * ...
                    obj.Cov_mat * Hmat');

                obj.Cov_mat = (eye(2 * obj.D_dim) - Kmat * Hmat) * obj.Cov_mat;
                Post_vec = [obj.State_vec; obj.Vel_vec] + Kmat * (Z_ass - ...
                    obj.State_vec);
                obj.State_vec = Post_vec(1 : obj.D_dim); % Z_ass
                obj.Vel_vec = Post_vec(obj.D_dim + 1 : end);
                % (obj.State_vec - obj.his_state(1, :)) / ...
                %     (obj.N_last * obj.T_ver)
                obj.his_state(end, :) = obj.State_vec';
            end
        end
    end
end

