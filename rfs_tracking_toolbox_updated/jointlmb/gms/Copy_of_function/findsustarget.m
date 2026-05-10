function [ x_sus,y_sus ] = findsustarget( multibeam_fm,cfar )
[Nx,Ny] = size(multibeam_fm);
x_sus = [];
y_sus = [];
[X,Y] = meshgrid(100:Nx-100,6:2*Ny-6);
[x_dectected,th]=cfar([multibeam_fm multibeam_fm],[X(:)';Y(:)']);
x_sus = [x_sus;X(x_dectected)];
y_sus = [y_sus;Y(x_dectected)];
end

