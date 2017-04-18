%Validation for measurements

addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Frenet')

%% generate a 3-D curve

%% add noise ? ------------------------------------------------------------
noise_strength=0;
noise_strength=1e-2;
rng(1); % set seed of random number generator

%% regularization weight --------------------------------------------------
weight=0; % default
weight=0.1; % regularizes torsion (but introduces small artifacts)

%% loop over examples 




%% specify arc-length sampling --------------------------------------------
ds=0.025; % space-step
slist=0:ds:10; % arc-length positions
ns=length(slist); % number of arc-length positions
lwin=50; % moving span

%% specify path curvature and torsion -------------------------------------

kappalist0=nan(size(slist));
taulist0=nan(size(slist));

        lambda=5; % wave-length
        kappalist0 =(slist/30); % path curvature
        taulist0(:)=cos(pi*slist/lambda); % path torsion        
            

%% integrate the Frenet-Seret equations to obtain the corresponding space curve
E0=eye(3); % identity matrix
% infinitesimal generators of rotations 
R1=[cross(E0(:,1),E0(:,1)) cross(E0(:,1),E0(:,2)) cross(E0(:,1),E0(:,3))]; % rotation around e1-axis
R2=[cross(E0(:,2),E0(:,1)) cross(E0(:,2),E0(:,2)) cross(E0(:,2),E0(:,3))]; % rotation aroudn e2-axis
R3=[cross(E0(:,3),E0(:,1)) cross(E0(:,3),E0(:,2)) cross(E0(:,3),E0(:,3))]; % rotation aroudn e3-axis

% allocate memory for co-moving Frenet frame and the path -----------------
E=nan(3,3,ns); % Frenet-frame in matrix notion as function of arc-length s
% E(:,3,is) ... tangent
% E(:,1,is) ... normal
% E(:,2,is) ... bi-normal
E(:,:,1)=E0; % set initial orientation
rr=nan(3,ns); % path as function of arc-length s
rr(:,1)=[0;0;0]; % set start point

% loop over all arc-length positions and do the integration ---------------
is=1;
for s=slist(1:end-1)
    kappa=kappalist0(is);
    tau=taulist0(is);
    R=kappa*R1+tau*R3; % local rotation of Frenet frame
    E(:,:,is+1)=E(:,:,is)*expm(R*ds); % propagate Frenet frame
    % E(:,:,is+1)=E(:,:,is)*(eye(3)+R*ds+0.5*R*R*ds^2); % same, but faster
    rr(:,is+1)=rr(:,is)+0.5*(E(:,3,is)+E(:,3,is+1))*ds; % propagate path along e3-vector
    is=is+1;
end; % s

%% add noise ==============================================================
%rr=rr+noise_strength*randn(size(rr));

%% Add radius

rad_curve= cos(pi/2*slist/lambda);
%rad_curve=rad_curve/max(max(rad_curve));    %make max r = 1;
%rad_curve=rad_curve./(3*kappalist);            % scale r so that the max is a will be at a thrid of R curvature















%% compute the Frenet frame ===============================================
[kappalist,taulist,ttlist,nnlist,bblist]=frenet_robust(rr,lwin,weight);
if nanmean(kappalist0.*kappalist')<0
    kappalist=-kappalist;
end;

% plot results ============================================================

figure, clf

% path
subplot(1,4,[1 2]), hold on
plot3(rr(1,:),rr(2,:), rr(3,:))
ind=1:10:ns;
quiver3(rr(1,ind),rr(2,ind), rr(3,ind),ttlist(1,ind),ttlist(2,ind),ttlist(3,ind),'g')
quiver3(rr(1,:),rr(2,:), rr(3,:),nnlist(1,:),nnlist(2,:),nnlist(3,:),'b')
quiver3(rr(1,:),rr(2,:), rr(3,:),bblist(1,:),bblist(2,:),bblist(3,:),'r')
daspect([1 1 1])
view([1 1 1])
title(sprintf('example no. %d',iexample))

% curvature
subplot(1,4,3), hold on
plot(slist,kappalist0,'r-','LineWidth',2) 
plot(slist,kappalist,'.-')
kappa_mean=nanmean(kappalist);
kappa_std=nanstd(kappalist);
% plot(slist,repmat(kappa_mean,ns,1),'g','LineWidth',2)
xlabel('s')
ylabel('\kappa')
ylim([min(min(kappalist0),kappa_mean-kappa_std)-0.1 max(max(kappalist0),kappa_mean-kappa_std)+0.1])
title('curvature \kappa')

% torsion
subplot(1,4,4), hold on
plot(slist,taulist0,'r-','LineWidth',2) 
plot(slist,taulist,'.-')
tau_mean=nanmean(taulist);
tau_std=nanstd(taulist);
xlabel('s')
ylabel('\tau')
ylim([min(min(taulist0),tau_mean-tau_std)-0.1 max(max(taulist0),tau_mean-tau_std)+0.1])
title('torsion \tau')

saveas(gcf,sprintf('example_%d.png',iexample))


