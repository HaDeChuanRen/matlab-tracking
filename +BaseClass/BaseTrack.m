classdef BaseTrack
% 基本的跟踪轨迹类，统一跟踪算法输出接口

    properties
        TrackID             % 跟踪ID int
        StartMoment         % 开始时间 double or datetime
        NLast               % 持续次数 int
        Dim                 % 跟踪信息维度 int
        TrackInfo           % 跟踪信息 TimeLast x Dim 矩阵
    end

    methods
        function obj = add(obj, newInfo)
            % 新增一行信息，持续时间加1
            % 输入: newInfo - 1 x Dim 的行向量，或 Dimx1 的列向量
            if nargin < 2
                error('新增TrackInfo需要提供新增信息行');
            end

            % 确保 newInfo 是行向量
            if iscolumn(newInfo)
                newInfo = newInfo';
            end

            % 检查维度
            if size(newInfo, 2) ~= obj.Dim
                error('新增信息的维度 (%d) 与 TrackInfo 的维度 (%d) 不匹配', ...
                    size(newInfo,2), obj.Dim);
            end

            % 添加新行
            obj.TrackInfo = [obj.TrackInfo; newInfo];
            obj.NLast = obj.NLast + 1;
        end

        function obj = update(obj, newInfo)
            % 修改最后一行的信息
            % 输入: newInfo - 1 x Dim 的行向量，或 Dimx1 的列向量
            if nargin < 2
                error('更新TrackInfo需要提供修改信息行');
            end

            % 确保 newInfo 是行向量
            if iscolumn(newInfo)
                newInfo = newInfo';
            end

            % 检查维度
            if size(newInfo, 2) ~= obj.Dim
                error('修改信息的维度 (%d) 与 TrackInfo 的维度 (%d) 不匹配', ...
                    size(newInfo,2), obj.Dim);
            end

            % 检查是否有最后一行
            if obj.NLast == 0
                error('当前轨迹没有信息行，无法修改');
            end

            % 修改最后一行
            obj.TrackInfo(end, :) = newInfo;
        end
    end

end

