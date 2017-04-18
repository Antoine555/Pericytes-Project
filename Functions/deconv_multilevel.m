% script for multilevel deconvolution

addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Multi level deconvolution\mltldeconvolution')

%% Generate data
% Original signal 
load I.mat
I = I(3:98,3:98,1:160);
edgeblur = fspecial('gaussian',5,4);
I = edgetaper(I,edgeblur);
%save ('I')

data.y_hat = fftn(I);

% PSF (9-by-9 uniform blur)
data.h_hat = zeros(size(I));
data.h_hat(1:9, 1:9, 1:9) = 1/729;
data.h_hat = circshift(data.h_hat, [-4 -4]);
data.h_hat = fftn(data.h_hat);

%set up debug
debug.xorig_hat = fftn(I);
debug.sigma = sqrt(2);

%% Algorithm parameters
% Wavelet family
% param.wavelet = 'Haar';
 param.wavelet = 'Daub2';
% param.wavelet = 'Sym4';
% param.wavelet = '97';

% Number of decomposition levels (for each dimension)
param.J = [3 3 3];

% Decimated or undecimated (redundant) transform
param.decimation = true;

% Random shifts
% param.shift = 'none';
param.shift = 'random';
% param.shift = 'cycle';

% Regularization parameter for the coarsest-scale scaling-function subband
param.lambdasf = 0;

% Regularization parameters for the wavelet subbands (level-dependent)
param.lambda = [0.04 0.04 0.04];

% Thresholding function
%param.threshold = 'Soft';
 param.threshold = 'Hard';

% Initial estimate
% param.iniest = 'zero';
 param.iniest = 'meas';
%param.iniest = 'tik';

% Number of iterations
param.K = 10;

% Type of multilevel cycles
%param.cycle = 'C'; % Coarse-to-fine
 param.cycle = 'V'; % V-cycles
% param.cycle = 'W'; % W-cycles
switch param.cycle
	case 'C'
		param.MG.eta1 = 0;
		param.MG.mu = 1;
		param.MG.eta2 = 1;
	case 'V'
		param.MG.eta1 = 1;
		param.MG.mu = 1;
		param.MG.eta2 = 1;
	case 'W'
		param.MG.eta1 = 1;
		param.MG.mu = 2;
		param.MG.eta2 = 1;
end

%% Run algorithm and compute SNR improvement
MLTLPrepare(data, param,debug);
result = MLTL(data, param,debug);
J=real(ifftn(result.x_hat));
%showimage(3, real(ifftn(result.x_hat)), 'Deconvolved');
%fprintf('SERG = %f\n', ComputeSERGain(debug.xorig_hat, data.y_hat, result.x_hat));
figure (1)
subplot(221);imshow(squeeze(I(:,50,:)),[0 255]);
axis equal
title('Original image');
subplot(222);imshow(squeeze(data.h_hat(:,:,round(size(data.h_hat,1)/2))),[]);
axis equal
title('xy psf');
subplot(223);imshow(squeeze(J(:,50,:)),[]);
axis equal
title('Deconvolved Image');
subplot(224);imshow(squeeze(data.h_hat(:,round(size(data.h_hat,1)/2),:)),[]);
axis equal
title('xz psf');
