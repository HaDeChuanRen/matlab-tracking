classdef CharaTarget
    %UNTITLED Summary of this class goes here
    % the amplitudes are not normalized, consider it

    properties
        ID_tar                  % target ID, int
        t_cnt                   % current time, double or datetime
        D_dim                   % dimension of point, int
        N_last                  % number of period which target lasts, int
        his_t                   % history of time, double
        T_ver                   % time duration of period, double

        num_points              % number of points, int
        point_set               % set of points, num_points * D_dim double
        amp_set                 % set of amplitudes, num_points * 1 double
        point_peak              % peak of all points, 1 * D_dim double
        Vel_vec                 % velocity of each dimension, 1 * D_dim double
        his_points              % history of target, N_last * 1 cell
        his_peak                % history of peak point, N_last * D_dim double
        P_s                     % probability of survive, double
    end

    methods
        function obj = CharaTarget(inCharaMea, inID, inT)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            if ~exist('inT', 'var'), inT = 1;
            elseif isempty(inT), inT = 1; end

            obj.ID_tar = inID;
            obj.D_dim = inCharaMea.D_dim;
            obj.num_points = inCharaMea.num_points;
            obj.point_set = inCharaMea.point_set;
            obj.amp_set = inCharaMea.amp_set;
            obj.Vel_vec = zeros(1, obj.D_dim);
            obj.T_ver = inT;
            obj.N_last = 1;
            obj.his_points = cell(1, 1);
            obj.his_points{1, 1} = obj.point_set;
            [~, peak_idx] = max(obj.amp_set);
            obj.point_peak = obj.point_set(peak_idx, :);
            obj.his_peak = obj.point_peak;
            obj.P_s = 0.3;
        end

        function obj = predict(obj, inT)
            if ~exist('inT', 'var'), inT = 1;
            elseif isempty(inT), inT = 1; end

            obj.T_ver = inT;
            obj.point_set = obj.point_set + obj.Vel_vec * inT;
            obj.P_s = obj.P_s - 0.1;
            obj.N_last = obj.N_last + 1;
            obj.his_points{end + 1, 1} = obj.point_set;

            [~, peak_idx] = max(obj.amp_set);
            obj.point_peak = obj.point_set(peak_idx, :);
            obj.his_peak = [obj.his_peak; obj.point_peak];
        end

        function obj = update(obj, Mea_ass)
            old_centroid = sum(obj.point_set .* obj.amp_set, 1)...
                ./ sum(obj.amp_set);
            mea_centroid = sum(Mea_ass.point_set .* Mea_ass.amp_set, 1)...
                ./ sum(Mea_ass.amp_set);
            obj.Vel_vec = (mea_centroid - old_centroid) / obj.T_ver;

            obj.num_points = Mea_ass.num_points;
            obj.point_set = Mea_ass.point_set;
            obj.amp_set = Mea_ass.amp_set;
            obj.his_points{obj.N_last, 1} = obj.point_set;
            obj.P_s = obj.P_s + 0.2;
            if (obj.P_s > 1), obj.P_s = 1; end

            [~, peak_idx] = max(obj.amp_set);
            obj.point_peak = obj.point_set(peak_idx, :);
            obj.his_peak(end, :) = obj.point_peak;
        end
    end
end