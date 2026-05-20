function [detSets, ampSets] = charaDets2D(shapeChar, detLoct, detSize)
    % charaDets_create 生成相应形状的观测列表
    arguments
        shapeChar (1, 1) char       % 形状字符
        detLoct (1, 2) double       % 观测中心
        detSize (1, 2) double       % 观测各维度尺寸
    end

    switch shapeChar
        case '.'
            % 点
            detSets = detLoct;
        case '|'
            % 竖条
            detSets = detLoct + [zeros(3, 1), detSize(2) * [-1/2; 0; 1/2]];
        case '-'
            % 横条
            detSets = detLoct + [detSize(1) * [-1/2; 0; 1/2], zeros(3, 1)];
        case '+'
            % 加号
            detSets = zeros(5, 2);
            detSets(1:3, :) = detLoct + [detSize(1) * [-1/2; 0; 1/2], ...
                zeros(3, 1)];
            detSets(4:5, :) = detLoct + [zeros(2, 1), detSize(2) * [-1/2; 1/2]];
        case 'd'
            % 菱形
            pointsNum = 4;
            thetaVec = (0 : pointsNum - 1)' * (2 * pi / pointsNum);
            detSets = detLoct + [detSize(1) * cos(thetaVec), ...
                detSize(2) * sin(thetaVec)];
        case 's'
            % 矩形
            detSets = detLoct + [detSize(1) * [1; -1; 1; -1], ...
                detSize(2) * [1; -1; -1; 1]];
        case '>'
            % 右三角
            pointsNum = 3;
            thetaVec = (0 : pointsNum - 1)' * (2 * pi / pointsNum);
            detSets = detLoct + [detSize(1) * cos(thetaVec), ...
                detSize(2) * sin(thetaVec)];
        case '<'
            % 左三角
            pointsNum = 3;
            thetaVec = (0 : pointsNum - 1)' * (2 * pi / pointsNum) + pi;
            detSets = detLoct + [detSize(1) * cos(thetaVec), ...
                detSize(2) * sin(thetaVec)];
        case '^'
            % 上三角
            pointsNum = 3;
            thetaVec = (0 : pointsNum - 1)' * (2 * pi / pointsNum) + pi/2;
            detSets = detLoct + [detSize(1) * cos(thetaVec), ...
                detSize(2) * sin(thetaVec)];
        case 'v'
            % 下三角
            pointsNum = 3;
            thetaVec = (0 : pointsNum - 1)' * (2 * pi / pointsNum) - pi/2;
            detSets = detLoct + [detSize(1) * cos(thetaVec), ...
                detSize(2) * sin(thetaVec)];
        case 'h'
            % 六边形
            pointsNum = 6;
            thetaVec = (0 : pointsNum - 1)' * (2 * pi / pointsNum);
            detSets = detLoct + [detSize(1) * cos(thetaVec), ...
                detSize(2) * sin(thetaVec)];
        case '*'
            % 六角星
            pointsNum = 6;
            thetaVec = (0 : pointsNum - 1)' * (2 * pi / pointsNum);
            detSets = detLoct + [detSize(1) * cos(thetaVec), ...
                detSize(2) * sin(thetaVec)];
            detSets = [detSets; detLoct];
        case '0'
            % 椭圆
            pointsNum = 8;
            thetaVec = (0 : pointsNum - 1)' * (2 * pi / pointsNum);
            detSets = detLoct + [detSize(1) * cos(thetaVec), ...
                detSize(2) * sin(thetaVec)];
        otherwise
    end
    pointsNum = size(detSets, 1);
    ampSets = ones(pointsNum, 1);

end