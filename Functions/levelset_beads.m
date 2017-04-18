% Level set script for the 3um beads in tissue and in water.
% If possible then a scripted wy of countunf the pixels along the central
% plane.

aspect_ratio=4.817734273;
wbeads=[763 93; 499 276; 441 538; 381 568; 114 695; ...
    944 715; 129 800; 826 958; 932 835];
window=30;



%% Stack non interpolated data
%First the tissue beads
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Capillaries work\Images\Bead FITC images\3 um perfused CTXa 2umStep imseq')
Stack_tbeads = [];     % Set empty matrix
%t=50/256;    % Set threshold
disp('Stacking started')
for i=0:1:8
    Stack_tbeads = cat(3,Stack_tbeads, imread(sprintf('C1-3 um perfused CTXa 2umStep0%d.tif',i)));
end

for i=10:1:45
    Stack_tbeads = cat(3,Stack_tbeads, imread(sprintf('C1-3 um perfused CTXa 2umStep%d.tif',i)));
end

%Second the water beads
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Capillaries work\Images\Calibration beads\240616 calibration bead imaging\3 um beads only take2 2umStep imseq')
Stack_wbeads = [];     % Set empty matrix
%t=50/256;    % Set threshold
disp('Stacking started')
for i=0:1:9
    Stack_wbeads = cat(3,Stack_wbeads, imread(sprintf('C1-3 um beads only take2 2umStep0%d.tif',i)));
end

for i=10:1:19
    Stack_wbeads = cat(3,Stack_wbeads, imread(sprintf('C1-3 um beads only take2 2umStep%d.tif',i)));
end
disp('Stacking finished and interpolation starting')


%% Interpolate the two stacks

Stack_tbeads_db=double(Stack_tbeads);
ny=size(Stack_tbeads_db,2);nx=size(Stack_tbeads_db,1);nz=size(Stack_tbeads_db,3)*4.817734273; %% desired output dimensions
[yq,xq,zq]=  ndgrid(linspace(1,size(Stack_tbeads_db,2),ny),...
    linspace(1,size(Stack_tbeads_db,1),nx),...
    linspace(1,size(Stack_tbeads_db,3),nz));
Stack_tbeads_int=interp3(Stack_tbeads_db,yq,xq,zq);


Stack_wbeads_db=double(Stack_wbeads);
ny=size(Stack_wbeads_db,2); nx=size(Stack_wbeads_db,1); nz=size(Stack_wbeads_db,3)*4.817734273; %% desired output dimensions
[yq,xq,zq]=  ndgrid(linspace(1,size(Stack_wbeads_db,2),ny),...
    linspace(1,size(Stack_wbeads_db,1),nx),...
    linspace(1,size(Stack_wbeads_db,3),nz));
Stack_wbeads_int=interp3(Stack_wbeads_db,yq,xq,zq);

disp('Interpolation finished and level set segmentation started.')


%% Level set the beads


addpath(genpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Level set toolbox\AOSLevelsetSegmentationToolboxM\data'));
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Level set toolbox\AOSLevelsetSegmentationToolboxM')

Stack_tbeads= padarray(Stack_tbeads(2:end-1,2:end-1,2:end-1),[1,1,1]);
Stack_wbeads= padarray(Stack_wbeads(2:end-1,2:end-1,2:end-1),[1,1,1]);

%{
for bead=1:1:size(tbeads,1)
    
    V_tbead=Stack_tbeads(...
        tbeads(bead,1)-(window+1)/2:tbeads(bead,1)+(window+1)/2,...
        tbeads(bead,2)-(window+1)/2:tbeads(bead,2)+(window+1)/2,:);
    
    % Initialise parameters
    smooth_weight = 2;
    image_weight = 1e-4;
    delta_t = 1;
    margin = 10;
    %phi = zeros(size(V)); % alternate initialisation
    %phi(margin:end-margin, margin:end-margin, margin:end-margin) = 1;
    phi=Stack_tbeads_int>130;
    phi = ac_reinit(phi-.5);
    
    % Run ChanVese level set on tissue beads volume
    for i = 1:2
        phi = ac_ChanVese_model(Stack_tbeads_int, phi, smooth_weight, image_weight, delta_t, 1);
        
        % for plotting:
        % if exist('h','var') && all(ishandle(h)), delete(h); end
        % iso = isosurface(phi);
        % h = patch(iso,'facecolor','w');  axis equal;  view(3);
        % set(gcf,'name', sprintf('#iters = %d',i));
        % view(3);
        % drawnow;
    end
    
    %Plot some slices to check result
    figure;
    slice = [30,35,40,45,50,55,60,65];
    for i = 1:8
        subplot(2,4,i); imshow(squeeze(Stack_tbeads(:,slice(i),:)),[]); hold on;
        c = contours(squeeze(phi(:,slice(i),:)),[0,0]);
        zy_plot_contours(c,'linewidth',2);
    end
    %find contour of deconvolved volume
    phi_seg=phi>0;
    
    figure;
    slice = [30,35,40,45,50,55,60,65];
    for i = 1:8
        subplot(2,4,i); imshow(squeeze(phi_seg(:,slice(i),:)),[]); hold on;
        c = contours(squeeze(phi(:,slice(i),:)),[0,0]);
        zy_plot_contours(c,'linewidth',2);
    end
    
    
    seg_imOut=(phi_seg);
    imOut1=Stack_tbeads;
    
    figure
    
    fv=isosurface(phi_seg,0);
    p = patch(fv);
    %isonormals(seg_stack3,p)
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    daspect([1,1,1])
    view(3); axis tight
    camlight
    lighting gouraud
    drawnow
    
    figure
    subplot(2,2,1); imagesc(Stack_tbeads(:,:,25))
    subplot(2,2,2); imagesc(squeeze(Stack_tbeads(:,50,:)))
    subplot(2,2,3); imagesc(phi_seg(:,:,25))
    subplot(2,2,4); imagesc(squeeze(phi_seg(50,:,:)))
    drawnow
    
    disp('level set of tissue beads finished')
end
%}

 % Do the same with the water beads
for bead=1:1:size(wbeads,1)
    
    V_wbead=Stack_wbeads_int(...
        wbeads(bead,1)-(window)/2:wbeads(bead,1)+(window)/2,...
        wbeads(bead,2)-(window)/2:wbeads(bead,2)+(window)/2,:); 
    
    % re-initialise parameters for level set of water beads
    smooth_weight = 2;
    image_weight = 1e-4;
    delta_t = 1;
    margin = 10;
    %phi = zeros(size(V)); % alternate initialisation
    %phi(margin:end-margin, margin:end-margin, margin:end-margin) = 1;
    phi2=V_wbead>130;
    phi2 = ac_reinit(phi2-.5);
    
    % run ChanVese level set on deconvolved volume
    for i = 1:4
        phi2 = ac_ChanVese_model(V_wbead, phi2, smooth_weight, image_weight, delta_t, 1);
        
        % for plotting:
        % if exist('h','var') && all(ishandle(h)), delete(h); end
        % iso = isosurface(phi);
        % h = patch(iso,'facecolor','w');  axis equal;  view(3);
        % set(gcf,'name', sprintf('#iters = %d',i));
        % view(3);
        % drawnow;
    end
    
    %Plot some slices to check result
    figure;
    imshow(squeeze(V_wbead(:,:,48)),[]); hold on;
        c2 = contours(squeeze(phi2(:,:,48)),[0,0]);
        zy_plot_contours(c2,'linewidth',2);
    
    %find contour of deconvolved volume
    phi_seg2=phi2>0;
    
    figure;
    imshow(squeeze(phi_seg2(:,:,48)),[]); hold on;
        c2 = contours(squeeze(phi2(:,:,48)),[0,0]);
        zy_plot_contours(c2,'linewidth',2);
       
    figure
    
    fv=isosurface(phi_seg2,0);
    p = patch(fv);
    %isonormals(seg_stack3,p)
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    daspect([1,1,1])
    view(3); axis tight
    camlight
    lighting gouraud
    drawnow
    
    figure
    subplot(2,2,1); imagesc(V_wbead(:,:,48))
    subplot(2,2,2); imagesc(squeeze(V_wbead(:,16,:)))
    subplot(2,2,3); imagesc(phi_seg2(:,:,48))
    subplot(2,2,4); imagesc(squeeze(phi_seg2(16,:,:)))
    drawnow
    disp('Segmentation finished.')
    
    %% Assessing the error
    bead_count=zeros(size(phi_seg2,3),1);
    for z_height=1:1:size(phi_seg2,3)
        count=0;
        for i=1:1:size(phi_seg2,1)
            for j=1:1:size(phi_seg2,2)
                if phi_seg2 (i,j,z_height) ==1
                    count=count+1;
                else
                end                              
            end
        end   
        bead_count(z_height)=count;
    end
    
    V_wbeads(bead)=struct('bead',bead,'stack',phi_seg2,'count',bead_count);
       
    
end


