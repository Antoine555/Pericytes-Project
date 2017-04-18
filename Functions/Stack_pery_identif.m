%Stacking Overlay for pericytes identification
clear all
close all
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Capillaries work\Images')
%% First Stack - magenta
% Load and import all perycites staining images into one 3D matrix
Stack = [];     % Set empty matrix
for i=1:1:9
    Stack = cat(3,Stack, imread(sprintf('26_z00%d_c001.tif',i)));
end

for i=10:1:42
   Stack = cat(3,Stack, imread(sprintf('26_z0%d_c001.tif',i)));
end

% define the SQUARE x-y window you would like to look at:
x_window= 51:274;   %pixels from top to bottom of images
y_window= 601:824;  %pixels from left to right of images
Stack_p= Stack(x_window,y_window,6:35);  % note I've set the slices to keep
% at 6th to 35th

S=size(Stack_p);
x=linspace (1,S(1),S(1));
y=linspace(1,S(2),S(2));
z=linspace (1,S(3),S(3));

%% Third Stack - green
% Load and import all vessel staining images into one 3D matrix
Stack3 = [];    % Set empty matrix
for i=1:1:9
    Stack3 = cat(3,Stack3, imread(sprintf('26_z00%d_c003.tif',i)));
end

for i=10:1:42
   Stack3 = cat(3,Stack3, imread(sprintf('26_z0%d_c003.tif',i)));
end

Stack3_p= Stack3(x_window,y_window,6:35);


%% Interpolation for magenta channel

Stack_pdb=double(Stack_p);
ny=size(Stack_pdb,2);nx=size(Stack_pdb,1);nz=size(Stack_p,3)*2;%4.817734273; %% desired output dimensions
[yq,xq,zq]=  ndgrid(linspace(1,size(Stack_p,2),ny),...
    linspace(1,size(Stack_p,1),nx),...
    linspace(1,size(Stack_p,3),nz));
imOut=interp3(y,x,z,Stack_pdb,yq,xq,zq);

imOut=permute(imOut, [2 1 3]);
S2=size(imOut);

%% Interpolation for green channel

Stack3_pdb=double(Stack3_p);
ny=size(Stack3_pdb,2);nx=size(Stack3_pdb,1);nz=size(Stack3_p,3)*2;%4.817734273; %% desired output dimensions
[yq,xq,zq]=  ndgrid(linspace(1,size(Stack3_p,2),ny),...
    linspace(1,size(Stack3_p,1),nx),...
    linspace(1,size(Stack3_p,3),nz));
imOut3=interp3(y,x,z,Stack3_pdb,yq,xq,zq);

imOut3=permute(imOut3, [2 1 3]);


%% Plotting
load MIP  % Maximum Intensity projection of magenta channel
x2=linspace (1,S2(1),S2(1));
y2=linspace(1,S2(2),S2(2));
z2=linspace (1,S2(3),S2(3));

figure
% Enter the parameters
% thresholds for isosurface; t1,t2 for perycites and t3 for vessels:
t1=70; t2=150; t3=130;  
% colors:
col1=[.7 .7 .8];        % Grey
col2=[1 .4 .4];         % Red
col3=[.3 1 .3];         % Green (for vessels)

% Plotting isosurfaces
hiso1 = patch(isosurface(x2,y2,z2,imOut,t1),'FaceColor',col1,...
    'EdgeColor','none');
axis equal;
view(3);
lighting phong;
isonormals(imOut,hiso1);
alpha(0.5);
set(gca,'DataAspectRatio',[1 1 1])
camlight;
drawnow
hold on;

hiso2 = patch(isosurface(x2,y2,z2,imOut,t2),'FaceColor',col2,...
    'EdgeColor','none');
isonormals(imOut,hiso2);
alpha(0.5);

hiso3 = patch(isosurface(x2,y2,z2,imOut3,t3),'FaceColor',col3,...
    'EdgeColor','none');
isonormals(imOut3,hiso3);
alpha(0.5);

imshow(imcomplement(MIP(x_window,y_window)));  % show MIP in plot