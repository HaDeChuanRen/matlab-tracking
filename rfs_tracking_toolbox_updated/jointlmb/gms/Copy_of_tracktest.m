clc;clear all;close all;fclose all;
addpath('.\Copy_of_function')
model= gen_model;
cfar = cfarInit(20,2);
histList= [];
meas.K = 0;
meas.Z = cell(101,1);
datadir = uigetdir('D:\1\识别专项数据\data\0402\20230402084303\bfout');
%%
datapath = dir([datadir,'\*.mat']);
if ~exist([datadir '\pic\'],'dir')
    mkdir([datadir '\pic\'])
end
for i = 50:150
    histList = clearList(histList);
    info = load([datadir '\' datapath(i).name]);
    [Nbeam,NFFT2] = size(info.multibeam_fm);
    signalType = checkSignalType( info.Tp,info.f0,info.type );
    switch signalType
        case 1
            copySignal = createCopySig(3e3,140,0.800,'hfm','up',NFFT2);
            startp = 1;
        case 2
            copySignal = createCopySig(3e3,140,0.400,'hfm','up',NFFT2);
            startp = 1500-816;
        case 3
            copySignal = createCopySig(3e3,140,0.400,'lfm','up',NFFT2);
            startp = 1500-816;
        case 4
            copySignal = createCopySig(2.73e3,140,0.400,'hfm','up',NFFT2);
            startp = 1500-816+512;
        case 5
            copySignal = createCopySig(3.28e3,140,0.400,'hfm','up',NFFT2);
            startp = 1500-816;
        otherwise
            error('未知信号')
    end
    disp(signalType)
    multibeam_mf = matchfilter(info.multibeam_fm.',copySignal.');
    multibeam_show = maxshow(multibeam_mf(startp:end,:),10,1);
    [susx,susy] = findsustarget(multibeam_show,cfar);
    addmeas = [];
    for w = 1:16
    if length(find(susy==i))>2
        susx(susy==w)=[];
        susy(susy==w)=[];
    end
    end
    if size(susx,1)>1
        d = pdist([susx,susy]);
        Z = linkage(d,'single');
        T = cluster(Z,'cutoff',10,'criterion','distance');
        for k = 1:max(T)
            tempx = susx(T==k);
            tempy = susy(T==k);
            tempy(tempy>16)=tempy(tempy>16)-16;
            index = tempx+size(multibeam_show,1)*(tempy-1);
            temp = multibeam_show(:);
            [~,idx] = max(temp(index,1));
            clusterX = tempx(idx);clusterY = tempy(idx);
            histList = updateList(clusterX,clusterY,histList);
            addmeas = [addmeas [100*clusterY;clusterX]];
        end
    elseif size(susx,1)==1
        histList = updateList(susx,susy,histList);
        
        addmeas = [addmeas [100*susy;susx]];
    end
    meas.K = meas.K+1;
    meas.Z{meas.K} = addmeas;
    [X,~,~] = run_filter_step(model,meas,meas.K);
    mainFig = figure(1);
    pcolor(abs([multibeam_show multibeam_show]));shading interp;colormap(jet);
    hold on
        plotList(histList)
    plot(X(1,:)/100,X(3,:),'ro','linewidth',2)
    hold off
    set(mainFig,'unit','centimeters','Position',[1 1 15 30])
    %     print(mainFig,'-dpng','-r0','-opengl',[datadir '\pic\' datapath(i).name '.png'])
    %     getpts();
end
%%
rmpath('.\Copy_of_function')