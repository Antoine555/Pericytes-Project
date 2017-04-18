% Testing of deconv_lucy on sphere

clear all; close all

%% Stacking
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Capillaries work\Images')
Stack3 = [];     % Set empty matrix
%t=50/256;    % Set threshold
for i=1:1:9
    Stack3 = cat(3,Stack3, imread(sprintf('26_z00%d_c003.tif',i)));
end

for i=10:1:42
    Stack3 = cat(3,Stack3, imread(sprintf('26_z0%d_c003.tif',i)));
end

%% Enter Perycites locations and ROI

peryxyzM=[600 600 20;183 199 29-4; 263 658 29-4];

window=150;
slice2D_size=20;
nber_slices= 25;
Excel_Data=zeros(nber_slices-4,size(peryxyzM,1)*5+1);
Excel_Data2=zeros(size(peryxyzM,1),2);

% Start the for loop to go through one perycite at a time
pery=1;
%store the current perycites location to use
peryxyz=peryxyzM(pery,:);

%% Take a new volume around the perycite
Stack3_p1= Stack3(peryxyz(1)-window/2:peryxyz(1)+window/2,peryxyz(2)-window/2:peryxyz(2)+window/2,5:38);
S=size(Stack3_p1);

x=linspace (1,S(1),S(1));
y=linspace(1,S(2),S(2));
z=linspace (1,S(3),S(3));

%% Use interpolation
%http://stackoverflow.com/questions/12520152/resizing-3d-matrix-image-in-matlab
Stack3_p1db=double(Stack3_p1);
ny=size(Stack3_p1db,2);nx=size(Stack3_p1db,1);nz=size(Stack3_p1,3)*4.817734273; %% desired output dimensions
[yq,xq,zq]=  ndgrid(linspace(1,size(Stack3_p1,2),ny),...
    linspace(1,size(Stack3_p1,1),nx),...
    linspace(1,size(Stack3_p1,3),nz));
imOut1=interp3(y,x,z,Stack3_p1db,yq,xq,zq);
S2=size(imOut1);

zi=peryxyz(3)*4.818-1;
peryxyzi=[peryxyz(1),peryxyz(2),zi]; %interpolate for pery z coord

x2=linspace (1,S2(1),S2(1));
y2=linspace(1,S2(2),S2(2));
z2=linspace (1,S2(3),S2(3));

I=imOut1;

%% Deconvolution

%Generate PSF
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Gaussian kernel')
%PSF=[];
%PSF = nonIsotropicGaussianPSF([1,1,4],3);
%[yq,xq,zq]=  ndgrid(linspace(1,size(PSF,2),31),...
%    linspace(1,size(PSF,1),31),...
%    linspace(1,size(PSF,3),61));
%PSF2=interp3(PSF,yq,xq,zq); %*255/0.04
%PSF2=PSF2(:,:,16:46);
PSF2=ones(30,30,30);
% Taper edges
edgeblur = fspecial('gaussian',5,4);
I = edgetaper(I,edgeblur);
I=I/255;

% Run R-L on I to obtain estimate of IO
clear SNR; clear ISNR; clear J;
%first deconvolution
J=deconvblind(I,PSF2,30);
%normalise J for next iteration
J=J/max(max(max(J)));

addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\3D median filter')

% Segment image
J_med=double(ordfilt3D(J,14));
J_med=J_med/max(max(max(J_med))); 
J_seg=J_med>0.2;

% Display results
figure (1)
subplot(221);imshow(squeeze(I(:,50,:)),[]);
axis equal
title('Raw blurred data');
subplot(222);imshow(squeeze(J(:,50,:)),[]);
axis equal
title('deconvolved image');
subplot(223);imshow(squeeze(J_med(:,50,:)),[]);
axis equal
title('Median filtered Image');
subplot(224);imshow(squeeze(J_seg(:,50,:)),[]);
axis equal
title('Segmented Image');

figure(3)
subplot(221);
fv=isosurface(x2,y2,z2,imOut1,52);
p = patch(fv);
%isonormals(seg_stack3,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud
drawnow

subplot(222);
fv=isosurface(x2,y2,z2,J,0.2);
p = patch(fv);
%isonormals(seg_stack3,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud
drawnow

subplot(223);
fv=isosurface(x2,y2,z2,J_seg,0);
p = patch(fv);
%isonormals(seg_stack3,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud
drawnow