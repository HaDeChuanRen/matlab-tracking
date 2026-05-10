clc;clear all;close all;fclose all;
addpath('.\function')
cfar = cfarInit(20,2);
histList= [];
datadir = uigetdir('D:\1\识别专项数据\7月data');
%%
datapath = dir([datadir,'\*.mat']);
if ~exist([datadir '\pic\'],'dir')
    mkdir([datadir '\pic\'])
end
for i = 1:length(datapath)
    histList = clearList(histList);
    info = load([datadir '\' datapath(i).name]);
    [Nbeam,NFFT2] = size(info.multibeam_fm);
    signalType = checkSignalType( info.Tp,info.f0,info.type );
    switch signalType
        case 1
            copySignal = createCopySig(2.73e3,140,0.400,'hfm','up',NFFT2);
            startp = 0;
        case 2
            copySignal = createCopySig(3.28e3,140,0.100,'hfm','up',NFFT2);
            startp = 1500-816;
        case 3
            copySignal = createCopySig(3e3,140,0.200,'hfm','up',NFFT2);
            startp = 1500-816;
        case 4
            copySignal = createCopySig(3e3,140,0.400,'hfm','up',NFFT2);
            startp = 1500-816+512;
        case 5
            copySignal = createCopySig(3e3,140,0.400,'hfm','up',NFFT2);
            startp = 1500-816;
        otherwise
            continue
    end
    disp(signalType)
    multibeam_mf = matchfilter(info.multibeam_fm.',copySignal.');
    multibeam_show = maxshow(multibeam_mf,10,1);
    [susx,susy] = findsustarget(multibeam_show,cfar);
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
        end
    elseif size(susx,1)==1
        histList = updateList(susx,susy,histList);
    end
    mainFig = figure(1);
    pcolor(abs([multibeam_show multibeam_show]));shading interp;colormap(jet);
    hold on
    plotList(histList)
    hold off
    set(mainFig,'unit','centimeters','Position',[1 1 15 30])
    print(mainFig,'-dpng','-r0','-opengl',[datadir '\pic\' datapath(i).name '.png'])
end
%%
rmpath('.\function')