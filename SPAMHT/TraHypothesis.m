classdef TraHypothesis
% Last update: 2025/2/28
properties
    Tar_set     % 目标集合
    Clu_set     % 杂波集合
    No_tol      % 目标标号总数
    D_dim       % 目标维度
    P_hyp       % 假设概率
    t_cnt       % 当前时刻
    t_his       % 时刻历史
    T_ver       % 追踪时间间隔

    Tra_his     % 历史轨迹
    Sigma_vec   % 目标测量方差
    P_Gate      % 门限概率
    L_ite       % 迭代次数
    Tar_ban     % 当前时刻消亡的目标
    Max_hyp     % 最大追踪假设数量
    Num_mote    % 杂波晋升所需有效关联周期数
end

methods
    function obj = TraHypothesis(inTarSet, inCluSet, inNo_tol, inPhyp, ...
        intcn, inTarhis, inSigma, inP_G, inL_ite, inTver, inMnum, inNmote)

        % default values of parameters
        if ~exist('inP_G', 'var'), inP_G = 1e-5;
        elseif isempty(inP_G), inP_G = 1e-5; end
        if ~exist('inL_ite', 'var'), inL_ite = 5;
        elseif isempty(inL_ite), inL_ite = 5; end
        if ~exist('inTver', 'var'), inTver = 1;
        elseif isempty(inTver), inTver = 1; end
        if ~exist('inMnum', 'var'), inMnum = 12;
        elseif isempty(inMnum), inMnum = 12; end
        if ~exist('inNmote', 'var'), inNmote = 20;
        elseif isempty(inNmote), inNmote = 20; end

        obj.Tar_set = inTarSet;
        obj.Clu_set = inCluSet;
        obj.No_tol = inNo_tol;
        obj.P_hyp = inPhyp;
        obj.t_cnt = intcn;
        obj.t_his = obj.t_cnt;
        obj.Tra_his = inTarhis;
        obj.P_Gate = inP_G;
        obj.L_ite = inL_ite;
        obj.Sigma_vec = inSigma;
        obj.D_dim = size(inSigma, 1);
        obj.T_ver = inTver;
        obj.Max_hyp = inMnum;
        obj.Num_mote = inNmote;
    end

    function [obj, All_Hyp] = SPA_track(obj, Mea_now, inMoment, inSigma)
        if exist('inMoment', 'var')
            obj.t_cnt = inMoment;
        else
            obj.t_cnt = obj.t_cnt + obj.T_ver;
        end
        obj.t_his = [obj.t_his; obj.t_cnt];
        if exist('inSigma', 'var'), obj.Sigma_vec = inSigma; end

        K_tar = length(obj.Tar_set);
        C_clu = length(obj.Clu_set);
        M_mea = size(Mea_now, 1);

        km_ass = [];
        K_show = [];
        All_Hyp = [];
        if (K_tar + C_clu) == 0
            for m_idx = 1 : M_mea
                clu_m = DimTarget(obj.No_tol, obj.t_cnt, Mea_now(m_idx, :)', ...
                    obj.Sigma_vec);
                obj.Clu_set = [obj.Clu_set; clu_m]; %#ok<*AGROW> 
                obj.No_tol = obj.No_tol + 1;
            end
            All_Hyp = obj;
        elseif (K_tar + C_clu) > 0
            [obj.Tar_set, Tar_pre] = obj.SetPredict(obj.Tar_set);
            [obj.Clu_set, Clu_pre] = obj.SetPredict(obj.Clu_set);

            if M_mea == 0
                All_Hyp = obj;
                return;
            end

            K_his = length(obj.Tra_his);
            Tar_dead = zeros(K_his, obj.D_dim);
            for kh_idx = 1 : K_his
                tk_end = obj.Tra_his(kh_idx).his_t(end);
                t_dead = obj.t_cnt - tk_end;
                for d_idx = 1 : obj.D_dim
                    Vel_kd = obj.Tra_his(kh_idx).Vel_vec(d_idx);
                    State_kd = obj.Tra_his(kh_idx).State_vec(d_idx);
                    Tar_dead(kh_idx, d_idx) = State_kd + Vel_kd * t_dead;
                end
            end

            [ProKM_a, ~] = obj.SPAassociate(Tar_pre, Mea_now);
            [km_ass, Pro_ass, Num_Hpy] = obj.Match_Create(ProKM_a);

            if isempty(km_ass)
                hyp_h = obj.Hyp_Create(Mea_now, [], obj.P_hyp, ProKM_a, ...
                    Tar_dead, Clu_pre);
                All_Hyp = [All_Hyp; hyp_h];
            else
                for h_idx = 1 : length(Pro_ass)
                    hyp_h = obj.Hyp_Create(Mea_now, km_ass(:, h_idx), ...
                        Pro_ass(h_idx) * obj.P_hyp, ProKM_a, Tar_dead, ...
                        Clu_pre);
                    All_Hyp = [All_Hyp; hyp_h];
                end
            end
            [~, Pmax_idx] = max(Pro_ass);
            obj = All_Hyp(Pmax_idx);
        end
        1;
    end

    function [NewTraHyp] = Hyp_Create(obj, Mea_now, km_ass, inPro, Promat_a, ...
        Tar_dead, clu_state)
        % late update: 2025/2/27
        % Detailed explanation goes here

        [New_set, Mea_remain] = obj.TarUpdate(Mea_now, km_ass, Promat_a);
        [ProCM_a, ~] = obj.SPAassociate(clu_state, Mea_remain);

        New_clu = obj.Clu_set;
        New_his = obj.Tra_his;
        New_tol = obj.No_tol;
        C_clu = size(clu_state, 1);

        Promo_c = [];
        ban_c = [];
        Ass_m = [];
        if ~isempty(ProCM_a)
            for c_idx = 1 : C_clu
                [Pro_ass, cm_ass] = max(ProCM_a(c_idx, :));
                if (cm_ass == size(Mea_remain, 1) + 1)
                    New_clu(c_idx) = New_clu(c_idx).hardUpdate([], 0);
                    if New_clu(c_idx).P_s < 1e-3
                        ban_c(end + 1) = c_idx;
                    end
                else
                    pd = 0.7;
                    New_clu(c_idx) = New_clu(c_idx).hardUpdate(Mea_remain(...
                        cm_ass, :)', Pro_ass, pd);
                    Ass_m(end + 1) = cm_ass;
                end
                if sum(New_clu(c_idx).his_Val) >= obj.Num_mote
                    Promo_c(end + 1) = c_idx;
                end
            end
        else
            for c_idx = 1 : C_clu
                New_clu(c_idx) = New_clu(c_idx).hardUpdate([], 0);
                if New_clu(c_idx).P_s < 1e-3
                    ban_c(end + 1) = c_idx;
                end
            end
        end

        % relife the dead targets by promoted targets
        Tar_Promote = zeros(length(Promo_c), obj.D_dim);
        for c_idx = 1 : length(Promo_c)
            c_Pro = Promo_c(c_idx);
            Tar_Promote(c_idx, :) = New_clu(c_Pro).State_vec';
        end
        ProCD_a = obj.SPAassociate(Tar_Promote, Tar_dead, obj.Sigma_vec);
        Tar_relife = [];
        if ~isempty(ProCD_a)
            for c_idx = 1 : length(Promo_c)
                c_Pro = Promo_c(c_idx);
                [Pro_cd, dead_idx] = max(ProCD_a(c_idx, :));
                if (Pro_cd > 0.5) && (dead_idx <= size(Tar_dead, 1))
                    New_clu(c_Pro).No_tar = obj.Tra_his(dead_idx).No_tar;
                    N_birth = find(obj.t_his == obj.Tra_his(dead_idx).his_t(1));
                    N_dead = find(obj.t_his == obj.Tra_his(dead_idx).his_t(end));
                    N_relife = find(obj.t_his == New_clu(c_Pro).his_t(1));
                    N_sleep = N_relife - N_dead - 1;
                    New_clu(c_Pro).his_t = obj.t_his(N_birth : end);
                    New_clu(c_Pro).N_last = length(New_clu(c_Pro).his_t);
                    if N_sleep > 0
                        his_coml = [];
                        for d_idx = 1 : obj.D_dim
                            his_comld = ...
                                linspace(obj.Tra_his(dead_idx).his_state(end,...
                                d_idx), New_clu(c_Pro).his_state(1, d_idx), ...
                                N_sleep)';
                            his_coml = [his_coml, his_comld];
                        end
                        New_clu(c_Pro).his_state = ...
                            [obj.Tra_his(dead_idx).his_state; his_coml; ...
                            New_clu(c_Pro).his_state];
                        New_clu(c_Pro).his_Ps = ...
                            [obj.Tra_his(dead_idx).his_Ps; zeros(N_sleep, 1); ...
                            New_clu(c_Pro).his_Ps];
                        New_clu(c_Pro).his_Val = ...
                            [obj.Tra_his(dead_idx).his_Val; zeros(N_sleep, 1); ...
                            New_clu(c_Pro).his_Val];
                    else
                        New_clu(c_Pro).his_state = ...
                            [obj.Tra_his(dead_idx).his_state(1 : end + N_sleep, ...
                            :); New_clu(c_Pro).his_state];
                        New_clu(c_Pro).his_Ps = ...
                            [obj.Tra_his(dead_idx).his_Ps(1 : end + N_sleep); ...
                            New_clu(c_Pro).his_Ps];
                        New_clu(c_Pro).his_Val = ...
                            [obj.Tra_his(dead_idx).his_Val(1 : end + N_sleep); ...
                            New_clu(c_Pro).his_Val];
                    end
                    Tar_relife = [Tar_relife; dead_idx];
                end
            end
        end
        New_his(Tar_relife) = [];


        % Promotion and banish the clutters
        New_set = [New_set; New_clu(Promo_c)];
        New_clu(union(Promo_c, ban_c)) = [];

        % create new clutters by measurements
        % Mea_new = setdiff(Mea_remain, Mea_remain(Ass_m));
        Mea_new = Mea_remain;
        Mea_new(Ass_m, :) = [];
        for m_idx = 1 : size(Mea_new, 1)
            Tar_m = DimTarget(New_tol, obj.t_cnt, Mea_new(m_idx, :)', ...
                obj.Sigma_vec);
            New_tol = New_tol + 1;
            New_clu = [New_clu; Tar_m];
        end

        % banish the undetected targets
        K_ban = [];
        K_tar = length(New_set);
        state_all = zeros(K_tar, 2 * obj.D_dim);
        for k_idx = 1 : K_tar
            if New_set(k_idx).P_s <= 1e-5
                K_ban(end + 1) = k_idx;
                continue;
            end
            state_k = [New_set(k_idx).State_vec', New_set(k_idx).Vel_vec'];
            [norm_min, min_idx] = min(sum(abs(state_all(1 : k_idx, :) - state_k), 2));
            if (norm_min == 0)
                if (New_set(k_idx).N_last > New_set(min_idx).N_last)
                    K_ban(end + 1) = min_idx;
                else
                    K_ban(end + 1) = k_idx;
                end
            else
                state_all(k_idx, :) = state_k;
            end
        end

        obj.Tar_ban = New_set(K_ban);
        for tar_idx = 1 : length(obj.Tar_ban)
            if(sum(obj.Tar_ban(tar_idx).his_Val) > obj.Num_mote)
                New_his = [New_his; obj.Tar_ban(tar_idx)];
            end
        end
        New_set(K_ban) = [];

        % banish the clutters
        C_clu = length(New_clu);
        C_ban = [];
        for c_idx = 1 : C_clu
            if(New_clu(c_idx).P_s < 1e-5)
                C_ban = [C_ban; c_idx];
            end
        end
        New_clu(C_ban) = [];

        NewTraHyp = obj;
        NewTraHyp.Tar_set = New_set;
        NewTraHyp.Clu_set = New_clu;
        NewTraHyp.Tra_his = New_his;
        NewTraHyp.No_tol = New_tol;
        NewTraHyp.P_hyp = inPro;
    end

    function [Promat_a, Promat_b] = SPAassociate(obj, tar_state, mea_now, inSig)
        if exist("inSig", 'var')
            Sig_vec = inSig;
        else
            Sig_vec = obj.Sigma_vec;
        end
        K_tar = size(tar_state, 1);
        M_mea = size(mea_now, 1);
        if isempty(tar_state) || isempty(mea_now)
            Promat_a = [];
            Promat_b = [];
            return;
        end
        Smat = zeros(K_tar, M_mea);
        for d_idx = 1 : obj.D_dim
            cost_d = abs(tar_state(:, d_idx) - mea_now(:, d_idx)');
            Smat = Smat + (cost_d .^ 2) / (Sig_vec(d_idx) ^ 2);
        end
        betamat = exp(-Smat);
        betamat(betamat < obj.P_Gate) = 0;
        udet_vec = obj.P_Gate * ones(K_tar, 1);
        betamat = [betamat, udet_vec];
        ximat = [ones(K_tar, M_mea); 0.01 * ones(1, M_mea)];

        % SPA algorithm iterate
        phimat_itel = betamat(:, 1 : M_mea) ./ (betamat(:, 1 + M_mea) + ...
            eps);
        vmat_itel = zeros(K_tar, M_mea);
        for l_iter = 1 : obj.L_ite
            vmat_itel = ximat(1 : K_tar, :) ./ (ximat(K_tar + 1, :) + ...
                sum(phimat_itel .* ximat(1 : K_tar, :), 1) - ...
                phimat_itel .* ximat(1 : K_tar, :));
            phimat_itel = betamat(:, 1 : M_mea) ./ (betamat(:, ...
                M_mea + 1) + sum(betamat(:, 1 : M_mea) .* vmat_itel, 2) ...
                - betamat(:,1 : M_mea) .* vmat_itel);
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
    end

    function [Match_mat, Match_pro, Match_Num] = Match_Create(obj, Promat_a, ...
        Pro_cut)
        Match_mat = [];
        Match_pro = 1;
        Match_Num = 1;

        if isempty(Promat_a)
            return;
        end

        if ~exist("Pro_cut", 'var'), Pro_cut = 1e-3;
        elseif isempty(Pro_cut), Pro_cut = 1e-3; end
        [K_tar, M_mea] = size(Promat_a);
        M_mea = M_mea - 1;

        Ass1_mat = Promat_a > Pro_cut;
        for k_idx = 1 : K_tar
            m_ass = find(Ass1_mat(k_idx, :));
            p_ass = Promat_a(k_idx, m_ass);
            m_ass(m_ass == M_mea + 1) = 0;
            Hnum_k = length(m_ass);
            Match_mat = repmat(Match_mat, [1 Hnum_k]);
            Match_pro = repmat(Match_pro, [1 Hnum_k]);
            km_k = zeros(1, Hnum_k * Match_Num);
            for h_idx = 1 : Hnum_k
                km_k((h_idx - 1) * Match_Num + (1 : Match_Num)) = m_ass(h_idx);
                Match_pro((h_idx - 1) * Match_Num + (1 : Match_Num)) = ...
                    Match_pro((h_idx - 1) * Match_Num + (1 : Match_Num)) * ...
                    p_ass(h_idx);
            end
            Match_mat = [Match_mat; km_k];

            prn_idx = find(Match_pro < Pro_cut);
            Match_mat(:, prn_idx) = [];
            Match_pro(prn_idx) = [];
            Match_pro = Match_pro / sum(Match_pro);
            Match_Num = size(Match_mat, 2);
        end

        if Match_Num == 0
            error("the number of Hypotheses is 0, please check it");
        end
        if (Match_Num > obj.Max_hyp)
            [~, Psort_idx] = sort(Match_pro, 'descend');
            Match_pro = Match_pro(Psort_idx(1 : obj.Max_hyp));
            Match_pro = Match_pro / sum(Match_pro);
            Match_mat = Match_mat(:, Psort_idx(1 : obj.Max_hyp));
        end
    end

    function [Dim_set, State_Pre] = SetPredict(obj, Dim_set, inSigma)
        if exist('inSigma', 'var'), obj.Sigma_vec = inSigma; end
        K_tar = length(Dim_set);
        State_Pre = zeros(K_tar, obj.D_dim);
        for k_idx = 1 : K_tar
            [Dim_set(k_idx), state_k] = Dim_set(k_idx).Predict(obj.t_cnt, ...
                obj.Sigma_vec);
            State_Pre(k_idx, :) = state_k';
        end
    end

    function [New_set, Mea_remain] = TarUpdate(obj, Mea_now, Match_vec, ...
            Promat_a, Pd_vec)
        New_set = [];
        km_del = Match_vec(Match_vec > 0);
        % Mea_remain = setdiff(Mea_now, Mea_now(km_del, :));
        Mea_remain = Mea_now;
        Mea_remain(km_del, :) = [];
        if isempty(obj.Tar_set)
            return; 
        end
        K_tar = length(obj.Tar_set);

        if ~exist("Pd_vec", 'var'), Pd_vec = 0.7 * ones(K_tar, 1);
        elseif isempty(Pd_vec), Pd_vec = 0.7 * ones(K_tar, 1); end

        for km_idx = 1 : K_tar
            m_ass = Match_vec(km_idx); 
            if (m_ass == 0)
                New_set = [New_set; obj.Tar_set(km_idx).hardUpdate([], 0)] ;
            else
                pro_ass = Promat_a(km_idx, m_ass);
                New_set = [New_set; obj.Tar_set(km_idx).hardUpdate(Mea_now(...
                    m_ass, :)', pro_ass,Pd_vec(km_idx))];
            end
        end
        
    end
end
end