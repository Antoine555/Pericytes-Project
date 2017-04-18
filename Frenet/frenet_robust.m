function [kappalist,taulist,ttlist,nnlist,bblist,scorelist]=frenet_robust(rr,lwin,weight)
%% computes robust estimates for path curvature and torsion using coalescing planes and circles
%
% INPUT: 
% rr ... 3xn-array of coordinates
% lwin ... span of moving windows: sets the spatial resolution for curvature and torsion
% weight ... a non-zero weight enfores continuity of Frenet frame [default=0, range=0..1]
%
% OUTPUT :
% kappalist ... nx1-array for curvature boundary values padded with NaNs
% taulist ... torsion
% ttlist ... 3xn array for tangent
% nnlist ... 3xn array for normal vector
% bblist ... 3xn array for binormal
% scorelist ... goodness of fit (0=perfect, >0.1 expect some problems, >0.2 failed)

% BM Friedrich, 08.09.2014

n=size(rr,2); % number of data points

% control moving windows --------------------------------------------------
nwin=n-lwin+1; % total number of moving windows

% control plot behavior ---------------------------------------------------
show_plot=true;
%show_plot=false;

% arc-length spacing of data points ---------------------------------------
ds_raw=sqrt(sum(gradient(rr,1,2).^2,1)); % arc-length spacing between control points
ds_smooth=smooth(ds_raw,round(0.4*lwin)); % mild smoothing increasing robustness in presence of noise
slist=cumsum(ds_smooth);

% useful functions --------------------------------------------------------
% normalize vector to unit length
normalize=@(vec) vec/norm(vec);

% alloc mem ---------------------------------------------------------------
ttlist=nan(size(rr)); % tangent
nnlist=nan(size(rr)); % normal
bblist=nan(size(rr)); % binormal
kappalist=nan(n,1); % signed path curvature
taulist=nan(n,1); % path torsion
scorelist=nan(n,1); % quality score for fitting coalescing plane

% loop over moving windows ================================================
for iwin=1:nwin  
    
    % select data points corresponding to current window ------------------
    rrwin=rr(:,iwin:lwin+iwin-1);
    ipos=iwin-1+round((1+lwin)/2); % corresponding index
    
    % fit straight line to estimate tangent -------------------------------
    fitmat=[ones(1,lwin);1:lwin]; % design matrix
    fit=fitmat'\rrwin'; % linear fit
    tt0=normalize(fit(2,:))'; % estimate for local tangent
    
    % fit coalescing plane ------------------------------------------------
    centroid=mean(rrwin,2);
    rrwin_centered=rrwin-repmat(centroid,1,lwin);
    % least-square fit of coalescing plane using singular-value-decomposition
    [U,S,~]=svd(rrwin_centered);
    % sort singular values
    [~,ind]=sort(diag(S),'descend');
    S=S(ind,ind);    
    U=U(:,ind);
    S=S/sum(diag(S)); % normalize singular values
    u1=U(:,1); % first unit vector spanning coalescing plane: agrees with tt0
    u2=U(:,2); % second unit vectors spanning coalescing plane: estimate for nn
    u3=U(:,3); % normal vector of coalescing plane
    
    % refine fit of coalescing plane: use previous fit as Bayesian prior --
    if (iwin>1) && (weight>1e-10)
        
        % cost function for fitting a plane with a Bayesian prior
        rotate_u2=@(phi) +cos(phi)*u2+sin(phi)*u3;
        rotate_u3=@(phi) -sin(phi)*u2+cos(phi)*u3;
        sse_fit_plane=@(phi,vec_prior,weight) ...
            ( cos(phi)*S(3,3) )^2 + ( sin(phi)*S(2,2) )^2+...
            -weight*vec_prior'*rotate_u3(phi);
        
        % start vector from SVD (without prior)
        phi_start=0;
                
        % specify Bayesian prior vector        
        vec_prior=bblist(:,ipos-1);

        % nonlinear fit of coalescing plane with Bayesian prior
        phi=fminsearch( @(phi) sse_fit_plane(phi,vec_prior,weight),phi_start);        
        
        % rotate u2 and u3 by angle phi around axis u1
        u2=rotate_u2(phi);
        u3=cross(u1,u2);
    end;
        
    % fit circle ----------------------------------------------------------
    param=CircleFitByTaubin([rrwin_centered'*u1 rrwin_centered'*u2]);
    xm=param(1); ym=param(2); % center of coalescing circle in (u1,u2)-plane
    RR=centroid+xm*u1+ym*u2; % center of coalescing circle in 3d-space
    r0=abs(param(3)); % radius of circle    

    % Frenet frame --------------------------------------------------------
    % phase angles corresponding to data points
    philist=sort(unwrap(atan2(rrwin_centered'*u2-ym,rrwin_centered'*u1-xm)));
    phim=angle(mean(exp(1i*philist))); % circular mean

    nn= u1*cos(phim)+u2*sin(phim); % normal
    tt=-u1*sin(phim)+u2*cos(phim); % tangent

    % tangent tt should point along the direction in which the path is transvered
    if tt0'*tt<0
        tt=-tt;
    end;
    
    % enforce continuity of Frenet frame ----------------------------------
    if iwin==1 
        % for the first window, specify trihedron orientation by hand     
        % - choose trihedron such that initial bb points along (+z)-direction
        if [0 0 1]*cross(tt,nn)<0            
            nn=-nn;
        end;
    else
        % all other windows: there should be no sign-flip wrt to previous window
        if nnlist(:,ipos-1)'*nn<0
            nn=-nn;
        end;
    end;
    bb=cross(tt,nn); % binormal
    
    % determine sign of signed curvature ----------------------------------
    if (centroid-RR)'*nn<0
        kappa=+1/r0;
    else
        kappa=-1/r0;
    end;
    
    % compute torsion from rotation of Frenet frame -----------------------
    if iwin==1
        tau=nan;
    else        
        dnn=+bblist(:,ipos-1)'*nn;        
        % dnn=-nnlist(:,ipos-1)'*bb; % gives the same result
        % simple formula, works fine for small ds
        % tau=dnn/(ds_smooth(ipos)*swin);
        % improved formula for finite ds
        psi=2*atan(dnn/2);
        tau=2*sin(psi/2)/ds_smooth(ipos);
    end;
    
    % keep a record -------------------------------------------------------    
    ttlist(:,ipos)=tt; % tangent
    nnlist(:,ipos)=nn; % normal
    bblist(:,ipos)=bb; % binormal 
    kappalist(ipos)=kappa; % signed curvature
    taulist(ipos)=tau; % torsion
    scorelist(ipos)=S(3,3)/S(2,2); % quality score for fitting coalescing plane

    % visual check --------------------------------------------------------
    if show_plot
        figure(1), clf, hold on
        plot3(rr(1,:),rr(2,:),rr(3,:),'color',0.5*[1 1 1])
        plot3(rrwin(1,:),rrwin(2,:),rrwin(3,:),'.-')
        plot3(centroid(1),centroid(2),centroid(3),'g.')
        plot3(centroid(1)+u1(1)*(xm+r0*cos(philist))+u2(1)*(ym+r0*sin(philist)),...
              centroid(2)+u1(2)*(xm+r0*cos(philist))+u2(2)*(ym+r0*sin(philist)),...
              centroid(3)+u1(3)*(xm+r0*cos(philist))+u2(3)*(ym+r0*sin(philist)),'r-');
        view(nn)  
        camup(bb)
        daspect([1 1 1])
    end;

end; % iwin
