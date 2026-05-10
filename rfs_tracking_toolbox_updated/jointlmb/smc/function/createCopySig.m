function [ cpyspec ] = createCopySig(f0,B,Tp,sigType,udflag,NFFT2)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
    fs = 20e3;
    fl=f0-B/2;
    fh = f0+B/2;
    m = B/Tp;
    D = 32;
    t=0:1/fs:Tp-1/fs;
    if strcmp('hfm',sigType)
        if strcmp('up',udflag)
            sourceFm=exp(-1i*2*pi*fl*(f0/m+Tp/2)*log(1-m*(t-Tp/2)/f0));
            bsourceFm=sourceFm(1:D:end).*exp(-1i*2*pi*f0*(1:D:length(sourceFm))/fs);
            cpyspec=conj(fft(bsourceFm,NFFT2));
        else
            sourceIFm=exp(1i*2*pi*fl*(f0/m+Tp/2)*log(1+m*(t-Tp/2)/f0));
            bsourceIFm=sourceIFm(1:D:end).*exp(-1i*2*pi*f0*(1:D:length(sourceIFm))/fs);   
            cpyspec=conj(fft(bsourceIFm,NFFT2));
        end
    else
        if strcmp('up',udflag)
            sourceFm=exp(1i*(2*pi*fl*t+pi*m*t.^2));
            bsourceFm=sourceFm(1:D:end).*exp(-1i*2*pi*f0*(1:D:length(sourceFm))/fs);   
            cpyspec=conj(fft(bsourceFm,NFFT2));
        else
            sourceIFm=exp(1i*(2*pi*fh*t-pi*m*t.^2));
            bsourceIFm=sourceIFm(1:D:end).*exp(-1i*2*pi*f0*(1:D:length(sourceIFm))/fs);   
            cpyspec=conj(fft(bsourceIFm,NFFT2));
        end
    end
end

