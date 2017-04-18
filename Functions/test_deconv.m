% Testing of deconv_lucy on sphere

clear all; close all
% Create a sphere inside wich elements have value 1 and 0 outside
n=linspace(-50,50,100);
[x,y,z]=ndgrid(n,n,n);
IO1=(sqrt((x+20).^2 + y.^2 + (z+20).^2)<=6);
IO2=(sqrt((x-17).^2 + y.^2 + (z-30).^2)<=3);
IO=IO1+IO2;

%Generate PSF
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Gaussian kernel')
PSF=[];
PSF = nonIsotropicGaussianPSF([1,1,4],3);
[yq,xq,zq]=  ndgrid(linspace(1,size(PSF,2),31),...
        linspace(1,size(PSF,1),31),...
        linspace(1,size(PSF,3),61));
PSF2=interp3(PSF,yq,xq,zq); %*255/0.04
PSF2=PSF2(:,:,16:46);

% Convolve IO with PSF 
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\3D convolution')
I_in = convolution3D_FFTdomain(IO,PSF2);
%IO_hat=fft(IO);
%PSF_hat=fft(PSF);
%I_hat=IO.*PSF;
%I=ifft(I_hat);

% Add Gaussian noise
I=I_in*255+5*randn(size(I_in));

% Taper edges
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
J_seg{1}=double(ordfilt3D(J{3},14)); %Use median filter for outliers
J_seg{2}=double(ordfilt3D(J{2},14));
J_seg{1}=J_seg{1}/max(max(max(J_seg{1})));  %Normalise again
J_seg{2}=J_seg{2}/max(max(max(J_seg{2})));
J_seg{1}=J_seg{1}>0.2;                      %Segment at 20%
J_seg{2}=J_seg{2}>0.2;

% Calculate SNRs with segmented images
SNR(2*(k-1)+1)=sum(sum(sum(J_seg{1}-IO)));
SNR(2*(k-1)+2)=sum(sum(sum(J_seg{2}-IO)));
noisy_psnr = 10 * log10(1/(sum(abs(IO(:) - I(:)).^2)/...
    size(I,1)/ size(I,2)));
PSNR(2*(k-1)+1) = 10 * log10(1/(sum(abs(IO(:) - J_seg{1}(:)).^2)/...
    size(I,1)/size(I,2)));
PSNR(2*(k-1)+2) = 10 * log10(1/(sum(abs(IO(:) - J_seg{2}(:)).^2)/...
    size(I,1)/size(I,2)));
ISNR(2*(k-1)+1) = PSNR(2*(k-1)+1) - noisy_psnr;
ISNR(2*(k-1)+2) = PSNR(2*(k-1)+2) - noisy_psnr;


%ISNR (2*(k-1)+1)= 10 * log10 ( sum( abs(IO(:) - I(:)).^2 ) ...
%                / sum( abs(IO(:) - J{3}(:)).^2 ) );
%ISNR (2*(k-1)+2)= 10 * log10 ( sum( abs(IO(:) - I(:)).^2 ) ...
%                / sum( abs(IO(:) - J{2}(:)).^2 ) );
%PSNR(2*(k-1)+1) = 10*log10(1/mean((IO(:)-J{3}(:)).^2));
%PSNR (2*(k-1)+2)= 10*log10(1/mean((IO(:)-J{2}(:)).^2));
k=k+1;
while k<=40
J=deconvlucy(J,PSF2,2);
%normalise J for next iteration
J{2}=J{2}/max(max(max(J{2})));
J{3}=J{3}/max(max(max(J{3})));

% Find the segmentation at this iteration
J_seg{1}=double(ordfilt3D(J{3},14)); %Use median filter for outliers
J_seg{2}=double(ordfilt3D(J{2},14));
J_seg{1}=J_seg{1}/max(max(max(J_seg{1})));  %Normalise again
J_seg{2}=J_seg{2}/max(max(max(J_seg{2})));
J_seg{1}=J_seg{1}>0.2;                      %Segment at 20%
J_seg{2}=J_seg{2}>0.2;

% Calculate SNRs with segmented images
SNR(2*(k-1)+1)=sum(sum(sum(J_seg{1}-IO)));
SNR(2*(k-1)+2)=sum(sum(sum(J_seg{2}-IO)));
noisy_psnr = 10 * log10(1/(sum(abs(IO(:) - I(:)).^2)/...
    size(I,1)/ size(I,2)));
PSNR(2*(k-1)+1) = 10 * log10(1/(sum(abs(IO(:) - J_seg{1}(:)).^2)/...
    size(I,1)/size(I,2)));
PSNR(2*(k-1)+2) = 10 * log10(1/(sum(abs(IO(:) - J_seg{2}(:)).^2)/...
    size(I,1)/size(I,2)));
ISNR(2*(k-1)+1) = PSNR(2*(k-1)+1) - noisy_psnr;
ISNR(2*(k-1)+2) = PSNR(2*(k-1)+2) - noisy_psnr;
k=k+1;
end

% Display results
figure (1)
subplot(221);imshow(squeeze(IO(:,50,:)),[]);
axis equal
title('Imaged object');
subplot(222);imshow(squeeze(I(:,50,:)),[]);
axis equal
title('Imaged object convolved with psf');
subplot(223);imshow(squeeze(J{2}(:,50,:)),[]);
axis equal
title('Deconvolved Image');
subplot(224);imshow(squeeze(PSF2(:,round(size(PSF2,1)/2),:)),[]);
axis equal
title('initial xz psf');

figure(2)
plot ([1:2*(k-1)],SNR)

figure (3)
plot ([1:2*(k-1)],ISNR)

figure (4)
plot ([1:2*(k-1)],PSNR)
