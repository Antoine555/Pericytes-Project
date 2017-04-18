% reconstruction BM4D script
clear all ;close all
%% Generate data
% Original signal
load I.mat
I=I(:,:,63:162);
edgeblur = fspecial('gaussian',5,4);
I = edgetaper(I,edgeblur);
%save ('I')

addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\BM4D')
% load constants
C = helper.constantsSparseTraj3D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modifiable parameters
%phantom    = C.SHEPPLOGAN;   % BRAINWEB SHEPPLOGAN
trajectory = C.RADIAL;       % RADIAL SPIRAL LOG_SPIRAL LIM_ANGLE SPHERICAL
%data       = C.COMPLEX;      % REAL COMPLEX

n          = 100;           % size of the 3d phantom (power of 2 <=256)
data_std   = 0;          % AWGN standard deviation in the initial observed data (%)
coverage   = 30.0;         % percentage of sampled pixels (%)
low_pass   = 9;            % number of retained phase coefficients (per dimension)
excursion  = 4;            % phase excursion (>=1)
min_norm   = eps;          % minimum normalized p-norm required to continue
pnorm      = 2;            % norm type used in early-stop condition

iter_nbr   = 10;          % max number of iterations
alpha      = 1.010;        % alpha noise-excitation
beta       = 5.0e2;        % beta noise-excitation

tol        = 1.0;          % tolerance error of coverage (%)
rot_deg    = 1*1e0;        % rotation degrees between consecutive trajectories
line_nbr   = 1;            % number of subsampling trajectories per slice
line_std   = 0.0;          % noise in subsampling trajectories (%)

lapse      = 10;           % lapse between saved slices during reconstruction
verbose    = C.IMAGE;        % NONE TEXT IMAGE
save_mat   = 0;            % save result to matlab .mat file

rec_type   = C.REC_3D;       % REC_2D REC_3D

do_wiener  = 1;            % perform BM4D Wiener filtering (1, 0)
profile    = 'np';         % BM4D parameter profile ('lc', 'np', 'mp')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       MODIFY BELOW THIS POINT ONLY IF YOU KNOW WHAT YOU ARE DOING       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% setting transform functions
% 3-D FFT
transform  = @fftn;
itransform = @ifftn;

% initial seeds
randn('seed',0);
rand('seed',0);




% getting sigma to be used in the first iteration
if data_std~=0
    sigma0 = data_std;
else
    sigma0 = 1000/eps;
end



% parameters' initialization

psnr_tilde      = zeros(1,iter_nbr);
psnr_hat        = zeros(1,iter_nbr);
ssim_tilde      = zeros(1,iter_nbr);
ssim_hat        = zeros(1,iter_nbr);
psnr_tilde_phi  = zeros(1,iter_nbr);
psnr_hat_phi    = zeros(1,iter_nbr);
ssim_tilde_phi  = zeros(1,iter_nbr);
ssim_hat_phi    = zeros(1,iter_nbr);
sw              = [1 1 1];

%buffer_cplx     = (y_hat_0.*exp(1i*phi_hat_0)) ./ sigma0^2;
%weight          = 1/sigma0^2;

%prev_y_tilde_k  = abs(buffer_cplx) ./ weight;

% initializing progression variable
idx = 2;

% Noise distribution
%distribution = 'Rice';
distribution = 'Gauss';

start = tic;
k = 1;
early_stop = 0;
y_hat_excite_k=I;
while k<=iter_nbr && ~early_stop
    % start counting
    iter = tic;
    
    if k>1
        
        % magnitude regularization
        y_hat_k = bm4d(y_hat_excite_k, distribution, 0, profile, do_wiener, 0);
        
        
        % buffer update
        %weight = weight + 1/sigma^2;
        %buffer_cplx = buffer_cplx + (y_hat_k.*exp(1i*phi_hat_k))/sigma^2;
        
        %y_tilde_k   = real(buffer_cplx./ weight);
        %phi_tilde_k = zeros(size(y_tilde_k));
        
        
        % early stop condition
        %early_stop      = (sum(abs(prev_y_tilde_k(:)-y_tilde_k(:)).^pnorm)).^(1/pnorm) / ...
        %    (numel(y_tilde_k)).^(1/pnorm) < min_norm;
        %prev_y_tilde_k  = y_tilde_k;
        
        
        y_hat_excite_k = y_hat_k;
        
        figure (1)
        subplot(121);imshow(squeeze(I(:,50,:)),[0 255]);
        axis equal
        subplot(122);imshow(squeeze(y_hat_k(:,50,:)),[]);
        axis equal
        title('Reconstructed Image');
        drawnow
       
    end
    
    % stop counting
    iter_time = toc(iter);
    
    k=k+1;
end
    
    
