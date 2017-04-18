% Richardson-Lucy deconvolution Script + denoising

% Set the random number generator back to its default settings for
% consistency in results.
rng default;
pause on
load I.mat
edgeblur = fspecial('gaussian',5,4);
I = edgetaper(I,edgeblur);
%save ('I')
%V = .0001;
%BlurredNoisy = imnoise(imfilter(I,PSF),'gaussian',0,V);
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\BM4D')
sigma             = 1;      % noise standard deviation given as percentage of the
                             % maximum intensity of the signal, must be in [0,100]
distribution      = 'Guass'; % noise distribution
                             %  'Gauss' --> Gaussian distribution
                             %  'Rice ' --> Rician Distribution
profile           = 'mp';    % BM4D parameter profile
                             %  'lc' --> low complexity
                             %  'np' --> normal profile
                             %  'mp' --> modified profile
                             % The modified profile is default in BM4D. For 
                             % details refer to the 2013 TIP paper.
do_wiener         = 0;       % Wiener filtering
                             %  1 --> enable Wiener filtering
                             %  0 --> disable Wiener filtering
verbose           = 1;       % verbose mode

estimate_sigma    = 0;       % enable sigma estimation

save_mat          = 0;       % save result to matlab .mat file
variable_noise    = 1;       % enable spatially varying noise
noise_factor      = 3;       % spatially varying noise range: [sigma, noise_factor*sigma]

%%
% check parameters
if sigma<=0
	error('Invalid "sigma" parameter: sigma must be greater than zero');
end
if noise_factor<=0
    error('Invalid "noise_factor" parameter: noise_factor must be greater than zero.');
end
estimate_sigma = estimate_sigma>0 || variable_noise>0;

% perform filtering
disp('Denoising started')
[y_est, sigma_est] = bm4d(I, distribution, (~estimate_sigma)*sigma, profile, do_wiener, verbose);

%Dispaly result
figure (1)
subplot(221);imshow(squeeze(I(:,50,:)),[0 255]);
axis equal
title('Original image');
subplot(222);imshow(squeeze(I(50,:,:)),[]);
axis equal
title('Original image 2');
subplot(223);imshow(squeeze(y_est(:,50,:)),[]);
axis equal
title('Deconvolved Image');
subplot(224);imshow(squeeze(y_est(50,:,:)),[]);
axis equal
title('Deconvolved Image 2')

pause

%%
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
INITPSF2=interp3(INITPSF,yq,xq,zq); %*255/0.04
INITPSF2=INITPSF2(:,:,16:46);

% Run the R-L deconvolution
J=deconvlucy(y_est,INITPSF2,25);

% Display results
figure (2)
subplot(221);imshow(squeeze(I(:,50,:)),[0 255]);
axis equal
title('Original image');
subplot(222);imshow(squeeze(INITPSF2(:,:,round(size(INITPSF2,1)/2))),[]);
axis equal
title('initial xy psf');
subplot(223);imshow(squeeze(J(:,50,:)),[]);
axis equal
title('Denoised and Deconvolved Image');
subplot(224);imshow(squeeze(INITPSF2(:,round(size(INITPSF2,1)/2),:)),[]);
axis equal
title('initial xz psf');

figure(3)
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