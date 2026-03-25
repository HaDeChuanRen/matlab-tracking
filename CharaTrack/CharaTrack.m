classdef CharaTrack < BaseClass.BaseTrack
    %UNTITLED Summary of this class goes here
    % the amplitudes are not normalized, consider it

    properties
        CntMoment               % 当前时刻, double or datetime
        MomentHis               % history of moments, double
        PeriodDuration            % time duration of period, double

        NumPoints               % number of points, int
        PointSet                % set of points, NumPoints * D_dim double
        AmpSet                  % set of amplitudes, NumPoints * 1 double
        PeakofPoints            % peak of all points, 1 * D_dim double
        Vels                    % velocity of each dimension, 1 * D_dim double
        PointsHis               % history of target, N_last * 1 cell
        PeakHis                 % history of peak point, N_last * D_dim double
        Ps                      % probability of survive, double
    end

    methods
        function obj = CharaTrack(inCharaMea, inID, inT)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            if ~exist('inT', 'var'), inT = 1;
            elseif isempty(inT), inT = 1; end

            obj.TrackID = inID;
            obj.Dim = inCharaMea.D_dim;
            obj.NumPoints = inCharaMea.num_points;
            obj.PointSet = inCharaMea.point_set;
            obj.AmpSet = inCharaMea.amp_set;
            obj.Vels = zeros(1, obj.Dim);
            obj.PeriodDuration = inT;
            obj.N_last = 1;
            obj.PointsHis = cell(1, 1);
            obj.PointsHis{1, 1} = obj.point_set;
            [~, peak_idx] = max(obj.amp_set);
            obj.PeakofPoints = obj.point_set(peak_idx, :);
            obj.PeakHis = obj.PeakofPoints;
            obj.Ps = 0.3;
        end

        function obj = predict(obj, inT)
            if ~exist('inT', 'var'), inT = 1;
            elseif isempty(inT), inT = 1; end

            obj.T_ver = inT;
            obj.point_set = obj.point_set + obj.Vels * inT;
            obj.Ps = obj.Ps - 0.1;
            obj.N_last = obj.N_last + 1;
            obj.PointsHis{end + 1, 1} = obj.point_set;

            [~, peak_idx] = max(obj.amp_set);
            obj.PeakofPoints = obj.point_set(peak_idx, :);
            obj.PeakHis = [obj.PeakHis; obj.PeakofPoints];
        end

        function obj = update(obj, Mea_ass)
            old_centroid = sum(obj.point_set .* obj.amp_set, 1)...
                ./ sum(obj.amp_set);
            mea_centroid = sum(Mea_ass.point_set .* Mea_ass.amp_set, 1)...
                ./ sum(Mea_ass.amp_set);
            obj.Vels = (mea_centroid - old_centroid) / obj.T_ver;

            obj.num_points = Mea_ass.num_points;
            obj.point_set = Mea_ass.point_set;
            obj.amp_set = Mea_ass.amp_set;
            obj.PointsHis{obj.N_last, 1} = obj.point_set;
            obj.Ps = obj.Ps + 0.2;
            if (obj.Ps > 1), obj.Ps = 1; end

            [~, peak_idx] = max(obj.amp_set);
            obj.PeakofPoints = obj.point_set(peak_idx, :);
            obj.PeakHis(end, :) = obj.PeakofPoints;
        end
    end
end