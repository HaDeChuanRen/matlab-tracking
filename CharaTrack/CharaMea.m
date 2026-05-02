classdef CharaMea
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Dim                           % dimension of point, int
        NumPoints                      % number of points, int
        PointSet                       % set of points, num_points * D_dim
        AmpSet                         % set of amplitudes, num_points * 1
    end

    methods
        function obj = CharaMea(inPointset, inAmpset)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            [obj.NumPoints, obj.Dim] = size(inPointset);
            obj.PointSet = inPointset;
            obj.AmpSet = inAmpset;
        end
    end
end