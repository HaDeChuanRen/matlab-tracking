classdef CharaMea
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        D_dim                           % dimension of point, int
        num_points                      % number of points, int
        point_set                       % set of points, num_points * D_dim
        amp_set                         % set of amplitudes, num_points * 1
    end

    methods
        function obj = CharaMea(in_pointset, in_ampset)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            [obj.num_points, obj.D_dim] = size(in_pointset);
            obj.point_set = in_pointset;
            obj.amp_set = in_ampset;
        end
    end
end