%Validation for measurements

addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Frenet')

%% generate a 3-D curve

%% add noise ? ------------------------------------------------------------
noise_strength=0;
noise_strength=1e-2;
rng(1); % set seed of random number generator


%% Generate straigth test tube

%% Create Test Tube
D=zeros(50,50,50);
k=1;
r=6+1*cos(k/6);
m=[35,35];
for k=16:50-15
    r=6+1*cos(k/6);
    for i=1:50
        for j=1:50
            if norm([i,j]-m)<=r
                D(i,j,k)=1;
            else
                D(i,j,k)=0;
            end
        end
    end
end

%theta=pi/4;
%[D_coord(:,1),D_coord(:,2),D_coord(:,3)]= ind2sub([50,50,50],find(D));
%D2_coord = D_coord*[0,cos(theta),sin(theta);...
%    0,-sin(theta),cos(theta);...
%    1,0,0];
%D2=zeros(50,50,50);
%for k=1:size(D2_coord)
%D2_ind=sub2ind([50,50,50],D2_coord(k,:));
%D2(D2_ind)=1;
%end

figure
x=linspace (1,size(D,1),size(D,1));
y=linspace (1,size(D,2),size(D,2));
z=linspace (1,size(D,3),size(D,3));

fv=isosurface(x,y,z,D,0.5);
p = patch(fv);
%isonormals(M,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight;
camlight
lighting gouraud
%[x,y,z]=ind2sub(D)
%plot3(D(1),D(2),D(3))



%% add noise ==============================================================
D=D+noise_strength*randn(size(D));


%% Add radius

%rad_curve= cos(pi/2*slist/lambda);
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


