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
        AmpsHis                 % history of peak point, N_last * D_dim double
        Ps                      % probability of survive, double
        Pd                      % probability of detection, double
    end

    methods
        function obj = CharaTrack(inCharaMea, inID, inTcnt, inTdur)
            % Construct an instance of CharaTrack class
            % Input:
            %   inCharaMea: input of characterized measurement, CharaMea
            %   inID: target ID, int
            %   inTcnt: current moment, double or datetime
            %   inTdur: duration of period, double or duration

            arguments
                inCharaMea (1,1) CharaMea
                inID (1,1) double
                inTcnt (1,1) = 1    % 可选：时间计数，默认1
                inTdur (1,:) = 1   % 可选：周期持续时间，默认1
            end

            % information of BaseTrack
            D_dim = inCharaMea.Dim;
            obj = obj@BaseClass.BaseTrack(inID, D_dim, inTcnt, zeros(1, D_dim));
            obj.TrackID = inID;
            obj.Dim = inCharaMea.Dim;
            obj.NLast = 1;
            obj.MomentInfo = inTcnt;

            % information of the character of targets
            obj.NumPoints = inCharaMea.NumPoints;
            obj.PointSet = inCharaMea.PointSet;
            obj.AmpSet = inCharaMea.AmpSet;
            obj.Vels = zeros(1, obj.Dim);
            obj.CntMoment = inTcnt;
            obj.PeriodDuration = inTdur;

            obj.PointsHis = cell(1, 1);
            obj.PointsHis{1, 1} = obj.PointSet;
            [~, peak_idx] = max(obj.AmpSet);
            obj.PeakofPoints = obj.PointSet(peak_idx, :);
            obj.AmpsHis{1, 1} = obj.AmpSet;

            % use the peak of points as the information of Track
            obj.TrackInfo = obj.PeakofPoints;
            obj.Ps = 0.3;
            obj.Pd = 0.98;
        end

        function obj = predict(obj, inTdur, weakRate)
            % predict the CharaTrack in the next moment
            if ~exist('inTdur', 'var') || isempty(inTdur), inTdur = 1; end
            if ~exist('weakRate', 'var') || isempty(weakRate), weakRate = 0.9; end

            obj.PeriodDuration = inTdur;
            obj.CntMoment = obj.CntMoment + inTdur;
            obj.PointSet = obj.PointSet + obj.Vels * inTdur;
            obj.Ps = obj.Ps * weakRate;  % 预测时生存概率衰减

            [~, peak_idx] = max(obj.AmpSet);
            obj.PeakofPoints = obj.PointSet(peak_idx, :);
            obj.PointsHis{end + 1, 1} = obj.PointSet;
            obj.AmpsHis{end + 1, 1} = obj.AmpSet;

            obj = obj.base_add(obj.PeakofPoints, obj.CntMoment);
        end

        function obj = update(obj, Mea_ass, assPro)
            old_centroid = sum(obj.PointsHis{obj.NLast - 1, 1} .* ...
                obj.AmpsHis{obj.NLast - 1, 1}, 1) ./ ...
                sum(obj.AmpsHis{obj.NLast - 1, 1});
            mea_centroid = sum(Mea_ass.PointSet .* Mea_ass.AmpSet, 1) ...
                ./ sum(Mea_ass.AmpSet);
            obj.Vels = (mea_centroid - old_centroid) / obj.PeriodDuration;

            obj.NumPoints = Mea_ass.NumPoints;
            obj.PointSet = Mea_ass.PointSet;
            obj.AmpSet = Mea_ass.AmpSet;
            obj.PointsHis{obj.NLast, 1} = obj.PointSet;
            obj.AmpsHis{obj.NLast, 1} = obj.AmpSet;
            % newPs = 1 - (1 - obj.Ps) * (1 - assPro) / (1 - obj.Pd * ...
            %     obj.Ps);
            obj.Ps = obj.Ps * (1 + assPro);
            if (obj.Ps >= 1), obj.Ps = 1 - 1e-6; end

            [~, peak_idx] = max(obj.AmpSet);
            obj.PeakofPoints = obj.PointSet(peak_idx, :);
            % obj.PeakHis(end, :) = obj.PeakofPoints;
            obj = obj.base_update(obj.PeakofPoints);
        end
    end
end