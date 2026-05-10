clc;clear all;close all
model= gen_model;
truth= gen_truth(model);
meas=  gen_meas(model,truth);
X= cell(meas.K,1);
N= zeros(meas.K,1);
L= cell(meas.K,1);
% for k= 1:meas.K
% [X{k},N(k),L{k}]=   run_filter_step(model,meas,k);
% end
est = run_filter(model,meas);
% est.X = X;est.N = N;est.L = L;
handles= plot_results(model,truth,meas,est);
% est = run_filter(model,meas);