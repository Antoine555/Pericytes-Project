% deconvolve with PSF stack before interpolation


%% Stacking PSF
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Capillaries work\Images\PSF images')
StackPSF = [];     % Set empty matrix
%t=50/256;    % Set threshold
for i=1:1:9
    StackPSF = cat(3,StackPSF, imread(sprintf('PSF BW0%d.tif',i)));
end

for i=10:1:42
    StackPSF = cat(3,StackPSF, imread(sprintf('PSF BW%d.tif',i)));
end

StackPSF=double(StackPSF(60:end-59,70:end-69,7:end-6));

I=Stack3_p1db;
PSF2=StackPSF;

%% Taper edges
edgeblur = fspecial('gaussian',5,4);
I = edgetaper(I,edgeblur);
I=I/255;

% Run R-L on I to obtain estimate of IO
clear SNR; clear ISNR; clear J;
k=1;
%first deconvolution
J=deconvlucy({I},PSF2,2);
%normalise J for next iteration
J{2}=J{2}/max(max(max(J{2})));
J{3}=J{3}/max(max(max(J{3})));

addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\3D median filter')
% Find the segmentation at this iteration
%J_seg{1}=double(ordfilt3D(J{3},14)); %Use median filter for outliers
%J_seg{2}=double(ordfilt3D(J{2},14));
%J_seg{1}=J_seg{1}/max(max(max(J_seg{1})));  %Normalise again
%J_seg{2}=J_seg{2}/max(max(max(J_seg{2})));
%J_seg{1}=J_seg{1}>0.2;                      %Segment at 20%
%J_seg{2}=J_seg{2}>0.2;


k=k+1;
while k<=40
    J=deconvlucy(J,PSF2,2);
    %normalise J for next iteration
    J{2}=J{2}/max(max(max(J{2})));
    J{3}=J{3}/max(max(max(J{3})));

    k=k+1;
end

% Only apply median filter 
J_med=double(ordfilt3D(J{2},14));
J_med=J_med/max(max(max(J_med))); 


%% Use interpolation
%http://stackoverflow.com/questions/12520152/resizing-3d-matrix-image-in-matlab
J_meddb=double(J_med);
ny=size(J_med,2);nx=size(J_med,1);nz=size(J_med,3)*4.817734273; %% desired output dimensions
[yq,xq,zq]=  ndgrid(linspace(1,size(J_med,2),ny),...
    linspace(1,size(J_med,1),nx),...
    linspace(1,size(J_med,3),nz));
imOutJ=interp3(J_meddb,yq,xq,zq);

imOutJ_seg=imOutJ>0.5;

% Display results
figure (1)
subplot(221);imshow(squeeze(I(:,50,:)),[]);
axis equal
title('Raw blurred data');
subplot(222);imshow(squeeze(J{2}(:,50,:)),[]);
axis equal
title('deconvolved image');
subplot(223);imshow(squeeze(J_med(:,50,:)),[]);
axis equal
title('Median filtered Image');
subplot(224);imshow(squeeze(imOutJ(:,50,:)),[]);
axis equal
title('Interpolated deconvolved filtered Image');

figure(2)
subplot(221);
fv=isosurface(imOut1,127);
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
fv=isosurface(J{2},0.5);
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
fv=isosurface(J_med,0.5);
p = patch(fv);
%isonormals(seg_stack3,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud
drawnow

subplot(224);
fv=isosurface(imOutJ,0.5);
p = patch(fv);
%isonormals(seg_stack3,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud
drawnow