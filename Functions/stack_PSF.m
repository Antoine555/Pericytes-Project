% Stack PSF
%% Stacking
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Capillaries work\Images\PSF images')
StackPSF = [];     % Set empty matrix
%t=50/256;    % Set threshold
for i=1:1:9
    StackPSF = cat(3,StackPSF, imread(sprintf('PSF BW0%d.tif',i)));
end

for i=10:1:42
    StackPSF = cat(3,StackPSF, imread(sprintf('PSF BW%d.tif',i)));
end

StackPSF=StackPSF((168/2+1-15):(168/2+1+15),(180/2+1-15):(180/2+1+15),...
    (42/2+1-15):(42/2+1+15));

%% Use interpolation
%http://stackoverflow.com/questions/12520152/resizing-3d-matrix-image-in-matlab
StackPSF_db=double(StackPSF);
ny=size(StackPSF,2);nx=size(StackPSF,1);nz=size(StackPSF,3)*4.817734273; %% desired output dimensions
[yq,xq,zq]=  ndgrid(linspace(1,size(StackPSF,2),ny),...
    linspace(1,size(StackPSF,1),nx),...
    linspace(1,size(StackPSF,3),nz));
imOutPSF=interp3(StackPSF_db,yq,xq,zq);

PSF2=imOutPSF(:,:,40:(end-40));
