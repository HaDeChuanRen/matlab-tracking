function [ multibeam_fm ] = matchfilter(multibeam,copysignal)
[Nx,Ny] = size(multibeam);
multibeam_fm = zeros(Nx,Ny);
for iy = 1:Ny
    multibeam_fm(:,iy)=abs(ifft(fft(multibeam(:,iy),Nx).*(copysignal)));
end
end

