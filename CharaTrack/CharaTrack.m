classdef CharaTrack < BaseClass.BaseTrack
    %UNTITLED Summary of this class goes here
    % the amplitudes are not normalized, consider it

    properties
        CntMoment               % current moment, double or datetime
        PeriodDuration          % time duration of period, double

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
        function obj = CharaTrack(inCharaMea, inID, inTcnt, inTdur)
            % Construct an instance of CharaTrack class
            % Input:
            %   inCharaMea: input of characterized measurement, CharaMea
            %   inID: target ID, int
            %   inTcnt: current moment, double or datetime
            %   inTdur: duration of period, double or duration

            if ~exist('inTcnt', 'var'), inTcnt = 1;
            elseif isempty(inTcnt), inTcnt = 1; end
            if ~exist('inTdur', 'var'), inTdur = 1;
            elseif isempty(inTdur), inTdur = 1; end

            % information of BaseTrack 
            obj.TrackID = inID;
            obj.Dim = inCharaMea.D_dim;
            obj.N_last = 1;
            obj.MomentInfo = inTcnt;

            obj.NumPoints = inCharaMea.num_points;
            obj.PointSet = inCharaMea.point_set;
            obj.AmpSet = inCharaMea.amp_set;
            obj.Vels = zeros(1, obj.Dim);
            obj.PeriodDuration = inTdur;
            
            obj.PointsHis = cell(1, 1);
            obj.PointsHis{1, 1} = obj.PointSet;
            [~, peak_idx] = max(obj.AmpSet);
            obj.PeakofPoints = obj.PointSet(peak_idx, :);
            obj.PeakHis = obj.PeakofPoints;
            obj.TrackInfo = obj.PeakofPoints;
            obj.Ps = 0.3;
        end

        function obj = predict(obj, inTdur)
            % predict the CharaTrack in the next moment
            if ~exist('inTdur', 'var') || isempty(inTdur), inTdur = 1; end

            obj.PeriodDuration = inTdur;
            obj.CntMoment = obj.CntMoment + inTdur;
            obj.PointSet = obj.PointSet + obj.Vels * inTdur;
            obj.Ps = obj.Ps - 0.1;

            [~, peak_idx] = max(obj.AmpSet);
            obj.PeakofPoints = obj.PointSet(peak_idx, :);
            obj.PeakHis = [obj.PeakHis; obj.PeakofPoints];
            obj.PointsHis{end + 1, 1} = obj.PointSet;

            obj.add(obj.PeakofPoints, obj.CntMoment);
        end

        function obj = update(obj, Mea_ass)
            old_centroid = sum(obj.PointSet .* obj.AmpSet, 1) ./ ...
                sum(obj.AmpSet);
            mea_centroid = sum(Mea_ass.point_set .* Mea_ass.amp_set, 1) ...
                ./ sum(Mea_ass.amp_set);
            obj.Vels = (mea_centroid - old_centroid) / obj.PeriodDuration;

            obj.NumPoints = Mea_ass.num_points;
            obj.PointSet = Mea_ass.point_set;
            obj.AmpSet = Mea_ass.amp_set;
            obj.PointsHis{obj.Nlast, 1} = obj.PointSet;
            obj.Ps = obj.Ps + 0.2;
            if (obj.Ps > 1), obj.Ps = 1; end

            [~, peak_idx] = max(obj.AmpSet);
            obj.PeakofPoints = obj.PointSet(peak_idx, :);
            % obj.PeakHis(end, :) = obj.PeakofPoints;
            obj.update(obj.PeakofPoints);
        end
    end
end