function [cfar] = cfarInit(NumTrainingCells,NumGuardCells)
% cfar = phased.CFARDetector('NumTrainingCells',NumTrainingCells,'NumGuardCells',NumGuardCells);
% exp_pfa = 1e-3;
% cfar.ThresholdFactor = 'Custom';
% cfar.CustomThresholdFactor=3;
% % cfar.ProbabilityFalseAlarm=exp_pfa;
% cfar.ThresholdOutputPort = true;
cfar = phased.CFARDetector2D('TrainingBandSize',[20,2], ...
    'GuardBandSize',[3,1],'ThresholdFactor','Custom','Method','GOCA', ...
    'CustomThresholdFactor',1.5,'ThresholdOutputPort',true);
end

