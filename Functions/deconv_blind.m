% Set the random number generator back to its default settings for
% consistency in results.
rng default;

I = imOut1(101:200,101:200,:);
edgeblur = fspecial('gaussian',5,4);
I = edgetaper(I,edgeblur);
%V = .0001;
%BlurredNoisy = imnoise(imfilter(I,PSF),'gaussian',0,V);

%Create a weight array to specify which pixels are included in processing
%WT = zeros(size(I));
%WT(5:end-4,5:end-4) = 1;
%INITPSF = ones(31,31,31);
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Gaussian kernel')
INITPSF=[];
INITPSF = nonIsotropicGaussianPSF([1,1,4],3);
[yq,xq,zq]=  ndgrid(linspace(1,size(INITPSF,2),31),...
        linspace(1,size(INITPSF,1),31),...
        linspace(1,size(INITPSF,3),61));
    INITPSF2=255/0.04*interp3(INITPSF,yq,xq,zq);
INITPSF2=INITPSF2(:,:,16:46);
%INITPSF =imgaussfilt3('gaussian',30,3);

% Perform blind deconvolution
[J,P] = deconvblind(I,INITPSF2,20,0.0001);

% Display results
figure (1)
subplot(221);imshow(squeeze(I(:,50,:)),[0 255]);
axis equal
title('Original image');
subplot(222);imshow(squeeze(P(:,:,round(size(P,1)/2))),[]);
axis equal
title('xy psf');
subplot(223);imshow(squeeze(J(:,50,:)),[]);
axis equal
title('Deconvolved Image');
subplot(224);imshow(squeeze(P(:,round(size(P,1)/2),:)),[]);
axis equal
title('xz psf');

figure(2)
subplot(121);
fv=isosurface(x2(101:200),y2(101:200),z2,imOut1(101:200,101:200,:),50);
p = patch(fv);
%isonormals(seg_stack3,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud
drawnow
subplot(122);
fv=isosurface(x2(1:100),y2(1:100),z2,J,50);
p = patch(fv);
%isonormals(seg_stack3,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud
drawnow