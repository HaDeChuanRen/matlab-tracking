classdef SPAMHTPlat
    %SPAMHTPLAT 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        hyp_all
        max_num
        d_dim
    end
    
    methods
        function obj = SPAMHTPlat(in_max, in_sigma, in_promote, in_time, ...
                in_pg, in_lite, in_tcir)
            %SPAMHTPLAT 构造此类的实例
            %   此处显示详细说明
            if ~exist('in_max', 'var'), in_max = 1;
            elseif isempty(in_max), in_max = 1; end
            if ~exist('in_sigma', 'var'), in_sigma = 1;
            elseif isempty(in_sigma), in_sigma = 1; end
            if ~exist('in_promote', 'var'), in_promote = 3;
            elseif isempty(in_promote), in_promote = 3; end
            if ~exist('in_time', 'var'), in_time = 0;
            elseif isempty(in_time), in_time = 0; end
            if ~exist('in_pg', 'var'), in_pg = 1e-5;
            elseif isempty(in_pg), in_pg = 1e-5; end
            if ~exist('in_lite', 'var'), in_lite = 5;
            elseif isempty(in_lite), in_lite = 5; end
            if ~exist('in_tcir', 'var'), in_tcir = 1;
            elseif isempty(in_tcir), in_tcir = 1; end
            

            obj.max_num = in_max;
            obj.d_dim = length(in_sigma);
            % obj.hyp_all = [];

            Tra_hyp = TraHypothesis([], [], 1, 1, in_time, [], in_sigma, ...
                in_pg, in_lite, in_tcir, obj.max_num, in_promote);
            obj.hyp_all = Tra_hyp;
        end
        
        function [obj, traj_all, target_all] = track(obj, mea_now)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明

            Hypset_snap = [];
            for h_idx = 1 : length(obj.hyp_all)
                [~, All_hyp] = obj.hyp_all(h_idx).SPA_track(mea_now);
                Hypset_snap = [Hypset_snap; All_hyp]; %#ok<AGROW> 
            end

            if length(Hypset_snap) > obj.max_num
                Pro_hyp = zeros(length(Hypset_snap), 1);
                for hp_idx = 1 : length(Hypset_snap)
                    Pro_hyp(hp_idx) = Hypset_snap(hp_idx).P_hyp;
                end
                [~, Psort_idx] = sort(Pro_hyp, 'descend');
                Pro_prn = Pro_hyp(Psort_idx(1 : Max_hyp));
                Pro_prn = Pro_prn / sum(Pro_prn);
                obj.hyp_all = Hypset_snap(Psort_idx(1 : obj.max_num));
                for h_idx = 1 : obj.max_num
                    obj.hyp_all(h_idx).P_hyp = Pro_prn(h_idx);
                end
            else
                obj.hyp_all = Hypset_snap;
            end

            target_all = [obj.hyp_all(1).Tar_set; obj.hyp_all(1).Tra_his];
            k_est = length(target_all);
            traj_all = cell(k_est, 2);
            for k_idx = 1 : k_est
                traj_all{k_idx, 1} = target_all(k_idx).his_state;
                traj_all{k_idx, 2} = target_all(k_idx).his_t;
            end
        end
    end
end

