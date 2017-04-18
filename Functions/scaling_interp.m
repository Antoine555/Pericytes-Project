clear all; close all

addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Capillaries work\Images')

Stack3 = [];     % Set empty matrix
%t=50/256;    % Set threshold
for i=1:1:9
    Stack3 = cat(3,Stack3, imread(sprintf('26_z00%d_c003.tif',i)));
end

for i=10:1:42
   Stack3 = cat(3,Stack3, imread(sprintf('26_z0%d_c003.tif',i)));
end

%portion of isosurface to be plotted and analysed
%pstart=601;             % starting pixel collumn at (pstart,pstart)
%pside=70;

Stack3_p= Stack3(1:180,1:180,1:30);

x=linspace (1,size(Stack3_p,1),size(Stack3_p,1));
y=linspace(1,size(Stack3_p,2),size(Stack3_p,2));
z=linspace (1,size(Stack3_p,3),size(Stack3_p,3));

% Plotting isosurfaces
t=30;    %threshold for isosurface
figure
subplot(1,2,1)
fv=isosurface(x,y,z,Stack3_p,t);
p = patch(fv); 
%isonormals(M,p)
p.FaceColor = 'green';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight 
lighting gouraud
drawnow

%% Use interpolation 
%http://stackoverflow.com/questions/12520152/resizing-3d-matrix-image-in-matlab
Stack3_pdb=double(Stack3_p);
ny=250;nx=250;nz=size(Stack3_p,3)*4.817734273; %% desired output dimensions
[yq,xq,zq]=  ndgrid(linspace(1,size(Stack3_p,1),ny),...
          linspace(1,size(Stack3_p,2),nx),...
          linspace(1,size(Stack3_p,3),nz));
imOut=interp3(Stack3_pdb,xq,yq,zq);

x2=linspace (1,size(imOut,1),size(imOut,1));
y2=linspace(1,size(imOut,2),size(imOut,2));
z2=linspace (1,size(imOut,3),size(imOut,3));

% Plotting isosurfaces
t=30;    %threshold for isosurface
subplot(1,2,2)
fv=isosurface(x2,y2,z2,imOut,t);
p = patch(fv); 
%isonormals(M,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight 
lighting gouraud
drawnow