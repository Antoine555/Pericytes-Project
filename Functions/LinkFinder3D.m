function [ output_args ] = LinkFinder3D( input_args )
%LINKFINDER3D: Function for automated selection and extraction of a 
% segmented link given a specified position, a stacked vessel volume (undecovolved
% & deconvolved) and a stacked red channel.
% It involves an interpolation for a specified window, a level set 3D 
% segmentation, a skeletonisation finally a pruning step.
% It outputs the coordinates of the centroids along the centreline as well
% as the gradient at these points. 


%Function for automated quantification of a vessel 
%metrics given a specified position, a stacked vessel volume (undecovolved
%& deconvolved) and a stacked red channel. 
%It involves an interpolation for a specified window, a level set 3D 
%segmentation, a skeletonisation and pruning, a 3D slicing function, a 2D
%level set segmentation. 
% It outputs the radii parameters, the enclosed area and the perimeter of
% the vessel on the vessel segment that is closest to the pericyte. An
% image projection is also generated.

% inputs: 
    % perixyzM : Matrix of coordinates for a set of points (nx3)
    % peri: pointer to the relevant row in perixyzM (scalar 1:n)
    % Stack_vessels: Imaged volume deconvolved but before interpolation
    %
    %    
    % window: 
    % slice2D_size: 
    % step_size: 
    % nber_slices_max:
    %
    %
    %
    
% Outputs: 
    % excluded: vector holding the pericytes that were too close to top/bot
    % 
    %
    %
    %
        
%% Initialisation
%store the current pericytes location to use
    perixyz=perixyzM(peri,:);
    perixyz(3)=round(perixyz(3)); %NB: this is the slice number!!!
    perixyz(1)=round(perixyz(1));
    perixyz(2)=round(perixyz(2));
    excluded=[];
    
    %% Take a volume around the pericyte
    disp('Checking pericyte position:')
    %check how close to the edge the pericyte is and act upon
    if  perixyz(3)<=5 || perixyz(3)>=(42-5)  ...
            || perixyz(1)<=20 || perixyz(1)>=(size(Stack_vessels,1)-20) ...
            || perixyz(2)<=20 || perixyz(2)>=(size(Stack_vessels,2)-20)
        % don't analyse this pericyte and go to the next
        excluded=cat(1,excluded,peri); % record which have been excluded
        disp(['pericyte ' num2str(peri) ' is too close to the top or'...
            ' bottom of sample. Going to the next pericyte.'])
        continue % go to next iteration in for loop
        % xe: coordinate of edges, x_peri_wind: peri within windowed Vol
    elseif perixyz(1)<= window/2
        xe1=1; x_peri_wind= perixyz(1);
        xe2=perixyz(1)+window/2;
        if perixyz(2)<= window/2
            ye1=1; y_peri_wind= perixyz(2);
            ye2=perixyz(2)+window/2;
            
        elseif perixyz(2)>=(size(Stack_vessels,2)-(window/2))
            ye2=size(Stack_vessels,2);
            y_peri_wind= window+1-size(Stack_vessels,2)+perixyz(2);
            ye1=perixyz(2)-window/2;
            
        else % perixyz(2) is far enough from the sides
            ye1=perixyz(2)-window/2;
            ye2=perixyz(2)+window/2;
            y_peri_wind= window/2+1;
            
        end
        
    elseif perixyz(1)>=(size(Stack_vessels,1)-(window/2))
        xe2=size(Stack_vessels,1);
        x_peri_wind= window+1-size(Stack_vessels,1)+perixyz(1);
        xe1=perixyz(1)-window/2;
        if perixyz(2)<= window/2
            ye1=1; y_peri_wind= perixyz(2);
            ye2=perixyz(2)+window/2;
            
        elseif perixyz(2)>=(size(Stack_vessels,2)-(window/2))
            ye2=size(Stack_vessels,2);
            y_peri_wind= window+1-size(Stack_vessels,2)+perixyz(2);
            ye1=perixyz(2)-window/2;
            
        else % perixyz(2) is far enough from the sides
            ye1=perixyz(2)-window/2;
            ye2=perixyz(2)+window/2;
            y_peri_wind= window/2+1;
            
        end
        
    else % perixyz(1) is far enough from the sides
        xe1=perixyz(1)-window/2;
        xe2=perixyz(1)+window/2;
        x_peri_wind= window/2+1;
        if perixyz(2)<= window/2
            ye1=1; y_peri_wind= perixyz(2);
            ye2=perixyz(2)+window/2;
            
        elseif perixyz(2)>=(size(Stack_vessels,2)-(window/2))
            ye2=size(Stack_vessels,2);
            y_peri_wind= window+1-size(Stack_vessels,2)+perixyz(2);
            ye1=perixyz(2)-window/2;
            
        else % perixyz(2) is far enough from the sides
            ye1=perixyz(2)-window/2;
            ye2=perixyz(2)+window/2;
            y_peri_wind= window/2+1;
            
        end
    end
    
    temp=[y_peri_wind, x_peri_wind];
    x_peri_wind=temp(1);
    y_peri_wind=temp(2);
    
    % take a portion of the non interpolated volume
    V0=Stack_vessels(xe1:xe2,ye1:ye2,1:end);
    Vref0=Vref00(xe1:xe2,ye1:ye2,1:end);
    S=size(V0);
    x=linspace (1,S(1),S(1));
    y=linspace(1,S(2),S(2));
    z=linspace (1,S(3),S(3));
    zi=(perixyz(3)-1)*4.818-1;  % This is the coordinate in the interpolated volume!!
    
    disp( ['Pericyte is valid and volume window has been selected.'...
        ' Interpolation of this volume started.'])
    
    %% interpolate this volume and call it V
    
    V0_db=double(V0);
    ny=size(V0,2);nx=size(V0,1);nz=size(V0,3)*4.817734273; %% desired output dimensions
    [yq,xq,zq]=  ndgrid(linspace(1,size(V0,2),ny),...
        linspace(1,size(V0,1),nx),...
        linspace(1,size(V0,3),nz));
    V=interp3(V0_db,yq,xq,zq);
    
    
    Vref0_db=double(Vref0);
    Vref=interp3(Vref0_db,yq,xq,zq);
    disp('Interpolation finished and level set segmentation started.')
    
    %% level set
    
    addpath(genpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Level set toolbox\AOSLevelsetSegmentationToolboxM\data'));
    addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Level set toolbox\AOSLevelsetSegmentationToolboxM')
    V= padarray(V(2:end-1,2:end-1,2:end-1),[1,1,1]);
    
    
    %%
    smooth_weight = 2;
    image_weight = 1e-4;
    delta_t = 1;
    margin = 10;
    %phi = zeros(size(V));
    %phi(margin:end-margin, margin:end-margin, margin:end-margin) = 1;
    phi=Vref>130;
    phi = ac_reinit(phi-.5);
    
    %%
    for i = 1:2
        phi = ac_ChanVese_model(V, phi, smooth_weight, image_weight, delta_t, 1);
        
        % if exist('h','var') && all(ishandle(h)), delete(h); end
        % iso = isosurface(phi);
        % h = patch(iso,'facecolor','w');  axis equal;  view(3);
        % set(gcf,'name', sprintf('#iters = %d',i));
        % view(3);
        % drawnow;
    end
    
    %%
    %{
        figure;
        slice = [30,35,40,45,50,55,60,65];
        for i = 1:8
            subplot(2,4,i); imshow(squeeze(V(:,slice(i),:)),[]); hold on;
            c = contours(squeeze(phi(:,slice(i),:)),[0,0]);
            zy_plot_contours(c,'linewidth',2);
        end
    %}
    phi_seg=phi>0;
    %{
        figure;
        slice = [30,35,40,45,50,55,60,65];
        for i = 1:8
            subplot(2,4,i); imshow(squeeze(phi_seg(:,slice(i),:)),[]); hold on;
            c = contours(squeeze(phi(:,slice(i),:)),[0,0]);
            zy_plot_contours(c,'linewidth',2);
        end
    %}
    seg_imOut=(phi_seg);
    imOut1=V;
    %{
        figure (1)
        
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
        
        figure(2)
        subplot(2,2,1); imagesc(Vref(:,:,25))
        subplot(2,2,2); imagesc(squeeze(Vref(:,50,:)))
        subplot(2,2,3); imagesc(phi_seg(:,:,25))
        subplot(2,2,4); imagesc(squeeze(phi_seg(50,:,:)))
        drawnow
    %}
    disp('Segmentation finished. Skeletonisation starting')
    %% Skeletonisation
    
    addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Skeleton3D')
    skel = Skeleton3D(seg_imOut);
    %{
        figure(3)
        col=[.7 .7 .8];
        hiso = patch(isosurface(permute(seg_imOut,[2 1 3]),0),'FaceColor',col,'EdgeColor','none');
        hiso2 = patch(isocaps(permute(seg_imOut,[2 1 3]),0),'FaceColor',col,'EdgeColor','none');
        axis equal;axis off;
        lighting phong;
        isonormals(permute(seg_imOut,[2 1 3]),hiso);
        alpha(0.5);
        set(gca,'DataAspectRatio',[1 1 1])
        camlight;
        hold on;
    %}
    w=size(skel,1);
    l=size(skel,2);
    h=size(skel,3);
    [x_skel,y_skel,z_skel]=ind2sub([w,l,h],find(skel(:)));
    %{
        plot3(x_skel,y_skel,z_skel,'.','Markersize',8,'MarkerFaceColor','b','Color','b');
        set(gcf,'Color','White');
        view(3)
        drawnow
    %}
    disp('Skeleton found. Pruning starting.')
    %% Cleaning and labelling the skeleton
    % Now use Skel2Graph on the skeleton
    addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Skel2graph3D V2')
    
    skel2=skel;
    
    w = size(skel2,1);
    l = size(skel2,2);
    h = size(skel2,3);
    
    % convert skeleton to graph structure
    [A,node,link] = Skel2Graph3D(skel2,17);
    
    % convert graph structure back to (cleaned) skeleton
    skel3 = Graph2Skel3D(node,link,w,l,h);
    
    i=1;
    % iteratively convert until there are no more 2-nodes left
    [A2,node2,link2] = Skel2Graph3D(skel3,17);
    while(min(cellfun('length',{node2.conn}))<3)
        skel3 = Graph2Skel3D(node2,link2,w,l,h);
        [A2,node2,link2] = Skel2Graph3D(skel3,17);
        i=i+1;
        if i>20
            break
        end
    end;
    %{
        % display result
        figure(4);
        hold on;
        for i=1:length(node2)
            x1 = node2(i).comx;
            y1 = node2(i).comy;
            z1 = node2(i).comz;
            xn(i)=x1;
            yn(i)=y1;
            zn(i)=z1;
            labels{i}={num2str(i)};
            for j=1:length(node2(i).links)    % draw all connections of each node
                if(node2(i).conn(j)<1)
                    col='b'; % branches are blue
                else
                    col='r'; % links are red
                end;
                
                % draw edges as lines using voxel positions
                for k=1:length(link2(node2(i).links(j)).point)-1
                    [x3,y3,z3]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k));
                    [x2,y2,z2]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k+1));
                    line([x3 x2],[y3 y2],[z3 z2],'Color',col,'LineWidth',3);
                end;
            end;
            
            % draw all nodes as yellow circles
            plot3(x1,y1,z1,'o','Markersize',9,...
                'MarkerFaceColor','y',...
                'Color','k');
        end;
        
        drawnow;
        view(-17,46);
        hold on
        col=[.7 .7 .8];
        hiso = patch(isosurface(permute(seg_imOut,[2 1 3]),0),'FaceColor',col,'EdgeColor','none');
        hiso2 = patch(isocaps(permute(seg_imOut,[2 1 3]),0),'FaceColor',col,'EdgeColor','none');
        axis equal; axis off;
        lighting phong;
        isonormals(permute(seg_imOut,[2 1 3]),hiso);
        alpha(0.5);
        view(3);
        %set(gca,'DataAspectRatio',[1 1 1])
        camlight;
        plot3 (x_peri_wind,y_peri_wind,zi,'MarkerSize',10,'MarkerFaceColor','m','MarkerEdgeColor','m','Marker','o')
    %}
    disp('Pruning finished. Finding closest link to pericyte.')
    %% Find distances of all centrepoints to pericyte
    
    indices=find(skel3(:));
    indices=squeeze(indices);
    distances=zeros(length(indices),1);
    for index=1:1:length(indices)
        
        %overallMinDistance = inf;
        refx=x_peri_wind;
        refy=y_peri_wind;
        refz=zi;
        [xtemp,ytemp,ztemp]=ind2sub(size(skel3),indices(index));
        distances(index) = sqrt((xtemp - refx).^2 + (ytemp - refy).^2+ (ztemp - refz).^2);
        
    end
    [minDistance, indexOfMin] = min(distances);
    %if minDistance < overallMinDistance
    % This boundary is the closer one.
    %end
    
    %% Find closest link
    closest_pt=indices(indexOfMin);
    
    % Now search for this point in the links
    for fields=1:1:length(link2);
        if sum(any(link2(fields).point==closest_pt))>=1;
            break
        elseif fields==length(link2);
            minDistance=100; % force the algo to stop analysing this pericyte
            % as the point identified is not on a 'link' but is on
            % the skeleton
        else
        end
    end
    
    % Extract this single link of the skeleton
    mylink=link2(fields).point;
    mynodes=[link2(fields).n1,link2(fields).n2];
    
    skel4=zeros(size(skel3));
    skel4(mylink)=1;   % Create the new matrix holding the single link
    S2=size(Vref);
    %figure(5)
    w=size(skel3,1);
    l=size(skel3,2);
    h=size(skel3,3);
    skel5=skel3-skel4;
    [x_skel5,y_skel5,z_skel5]=ind2sub([w,l,h],find(skel5(:)));
    %{
        plot3(x_skel5,y_skel5,z_skel5,'.','Markersize',8,'MarkerFaceColor','b','Color','b');
        axis ([0 S2(1) 0 S2(2) 0 S2(3)]);
        drawnow
        hold on
    %}
    w=size(skel4,1);
    l=size(skel4,2);
    h=size(skel4,3);
    [x_skel4,y_skel4,z_skel4]=ind2sub([w,l,h],find(skel4(:)));
    %{
        plot3(x_skel4,y_skel4,z_skel4,'.','Markersize',8,'MarkerFaceColor','r','Color','r');
        view(3)
        drawnow
        hold on
        axis tight; axis off;
        daspect([1,1,1])
        view(3);
        plot3 (window/2,window/2,zi,'MarkerSize',7,'MarkerFaceColor','m','MarkerEdgeColor','m','Marker','o')
    %}
    %% Continue the analyisis if the link found is correct
    %skip this pericyte if the correct link could not be computed
    if minDistance>=20 % link is too far from perictyte
        excluded=cat(1,excluded,peri);
        disp(['Pericyte ' num2str(peri) ' excluded as too far from'...
            ' nearest link. /n' 'Check z coordinate or increase window size. /n'...
            'Press any key to go to next peri.'])
        pause
        continue
    else
    end
    
    %% Use spline interpolation to smooth the centreline
    addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\SmoothN');
    %figure (6)
    mylink_coord=zeros(size(mylink,2),3);
    for ind=1:1:size(mylink,2)
        [mylink_coord(ind,1),mylink_coord(ind,2),mylink_coord(ind,3)]=...
            ind2sub(size(skel4),mylink(ind));
    end
    mylink_x=mylink_coord(:,1);
    mylink_y=mylink_coord(:,2);
    mylink_z=mylink_coord(:,3);
    
    t = cumsum([0;sqrt(diff(mylink_x).^2 + diff(mylink_y).^2 + diff(mylink_z).^2)]);
    tt = t(1):2:t(end);
    u = smoothn({mylink_x,mylink_y,mylink_z},100);
    %{
            plot3 (mylink_x,mylink_y,mylink_z,'.',u{1},u{2},u{3},'r')
            axis tight; axis off;
            daspect([1,1,1])
            view(3)
            box on
            pause
    %}
    disp('Valid link found and smoothed.')
    %% interpolate for equidistant points along the smooth centreline and diff
    %http://fr.mathworks.com/matlabcentral/fileexchange/34874-interparc
    addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\interparc')
    
    F=[u{1},u{2},u{3}]; %Matrix of coordinates for the smooth spline
    %interpolate to find arclength of link
    [F2,dF2dt] = interparc(3*size(mylink_x,1),mylink_x,mylink_y,mylink_z);
    mylink_al=norm(F2(3,:)-F2(2,:))*3*size(mylink_x,1);
    
    %interpolate to set slicing points
    nber_slices= round(mylink_al/step_size);     % number of equidistant slices to be taken
    [F3,dF3dt] = interparc(nber_slices,F(:,1),F(:,2),F(:,3)); % also returns derivatives, and a function handle
    
    
    % check if there are nodes on either side of the link
    mynodes2_coord=[0 0 0];
    if mynodes(2)<=0            % if not use the last point on the link
        [mynodes2_coord(1),mynodes2_coord(2),mynodes2_coord(3)] = ind2sub(size(skel),link2(fields).point(end));
        mylink_l=norm(...
            [node2(mynodes(1)).comx,node2(mynodes(1)).comy,node2(mynodes(1)).comz] ...
            -[mynodes2_coord(1),mynodes2_coord(2),mynodes2_coord(3)]);
    else
        mylink_l=norm(...
            [node2(mynodes(1)).comx,node2(mynodes(1)).comy,node2(mynodes(1)).comz] ...
            -[node2(mynodes(2)).comx,node2(mynodes(2)).comy,node2(mynodes(2)).comz]);
    end
    
    %% Find the closest point to pericyte on new smooth link
    distances=zeros(size(F3,1),1);
    for index=1:1:size(F3,1)
        
        %overallMinDistance = inf;
        refx=x_peri_wind;
        refy=y_peri_wind;
        refz=zi;
        distances(index) = sqrt((F3(index,1) - refx).^2 +...
            (F3(index,2) - refy).^2+ (F3(index,3) - refz).^2);
        
    end
    [minDistance, indexOfMin] = min(distances);
    
    if indexOfMin<=5 || size(F3,1)-indexOfMin<=6
        disp(['Pericyte ' num2str(peri) ' is too close to an endpoint'...
            ' (probably at a bifurcation. This script cannot yet'...
            ' analyse these. Please close all open windows then' ...
            ' press any key to go to the next pericyte.'])
        pause
        continue
    else
    end
    disp('Starting slicing at centrepoints.')
    
    %% Extract Slices
    
    %figure (7)
    % store value for the latest skeleton
    c_skel4=[x_skel4,y_skel4,z_skel4];
    %pause on
    map=[0 0 0; 0 0.2 0;0 0.4 0;0 0.6 0; 0 0.8 0;0 1 0];
    
    %Extract slices and plot them onto isosurface
    addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Extract Slice')
    Stored_radii_para=[];%set matrix for storing mean max min radii at each centrepoint
    d=[0 0 0]; %initialise direction vector
    %slice2D_size= 20; %choose parametric size of the slices taken
    
    slice2D_temp=zeros(slice2D_size*2+1,slice2D_size*2+1,size(F3,1));
    sliceInd_temp=zeros(slice2D_size*2+1,slice2D_size*2+1,size(F3,1));
    subX_temp=zeros(slice2D_size*2+1,slice2D_size*2+1,size(F3,1));
    subY_temp=zeros(slice2D_size*2+1,slice2D_size*2+1,size(F3,1));
    subZ_temp=zeros(slice2D_size*2+1,slice2D_size*2+1,size(F3,1));
    %c_2D=zeros(1,size(F3,1)-1); %centre of image is fixed for now
    
    for row=1:size(F3,1) %for all centrepoints
        % Find direction vector (more than one method)
        % d=F3(row,:)-F3(row-1,:); % note cannot start at row one
        % d=d*5;
        
        % Find direction vector using derivative at point
        d=dF3dt(row,:);
        d=d/max(d)*50;
        
        % extract slices if using the point on smooth spline as origin:
        [slice2D_temp(:,:,row),sliceInd_temp(:,:,row),...
            subX_temp(:,:,row),subY_temp(:,:,row),subZ_temp(:,:,row)]...
            =extractSlice(V,...
            F3(row,1),F3(row,2),F3(row,3),d(1),d(2),d(3),slice2D_size);
    end
    slice2D_temp(isnan(slice2D_temp)) = 0 ;
    
    for row=1:size(F3,1)
        %Plot slices as we circulate through the centrepoints
        surf(subX_temp(:,:,row),subY_temp(:,:,row),subZ_temp(:,:,row),...
            slice2D_temp(:,:,row),'FaceColor','texturemap','EdgeColor','none');
        axis([1 S2(1) 1 S2(2) 1 S2(3)]); colormap(map);
        view(-17,46);
        drawnow
        hold on
    end
    %{
    pause
    col=[.7 .7 .8];
    hiso = patch(isosurface(permute(seg_imOut,[2 1 3]),0),'FaceColor',col,'EdgeColor','none');
    hiso2 = patch(isocaps(permute(seg_imOut,[2 1 3]),0),'FaceColor',col,'EdgeColor','none');
    axis equal;%axis off;
    lighting phong;
    isonormals(permute(seg_imOut,[2 1 3]),hiso);
    alpha(0.5);
    %set(gca,'DataAspectRatio',[1 1 1])
    camlight;
    pause
    %}
    % Plot all the 2D slices
    figure
    for row=1:size(F3,1)
        subplot(round((size(F3,1))/8+1),8,row)
        imagesc(slice2D_temp(:,:,row))
    end
    
    pause(0.5)
    
    disp(['Slices found, now recalculating the contours' ...
        ' and centoids in 2D'])
    
    %% find the contour using hybrid level set
    
    addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Level set toolbox\AOSLevelsetSegmentationToolboxM')
    seg_slice2D_temp=zeros(slice2D_size*2+1,slice2D_size*2+1,size(F3,1)-4);
    contour_ptM=[];
    phi=[];
    contour_pt_s=struct('cross_section',0,'coord',[]);
    fprintf([ '# of slices (-4): ' num2str(size(F3,1)-4) ...
        '. Contours found:  '])
    for row=3:1:size(F3,1)-2
        if (row-3)<=10
            fprintf('\b%d',(row-3)); % update counter display
        else
            fprintf('\b\b%d',(row-3)); % update counter display
        end
        slice2D_temp(:,:,row)=slice2D_temp(:,:,row)./max(max(slice2D_temp(:,:,row)))*255;
        %I =  ordfilt3D(slice2D_temp(:,:,row),14);
        I =  slice2D_temp(:,:,row);
        %phi = ac_SDF_2D('circle', size(I), c_2D(15,:),slice2D_size) ;
        phi = ac_SDF_2D('circle', size(I), [round(size(I,1)/2),round(size(I)/2)] ,2);
        
        g = ac_gradient_map(I, 5);
        mu = 0.1*(slice2D_temp(16,16,row));
        propagation_weight = 0.6; GAC_weight = 5000;
        delta_t = 1; n_ite = 100; show_result = 0; n_iters=100;
        phi = ac_hybrid_model(I-mu, phi, propagation_weight, GAC_weight, g, ...
            delta_t, n_iters, show_result);
        %close all
        seg_slice2D_temp(:,:,row)=phi>=0;
        
        c = contours(phi,[0,0]);
        contour_pt=[];
        idx = 1;
        while idx < size(c,2)
            n = c(2,idx);
            contour_pt = [contour_pt, c(:,idx+1:idx+n)];
            idx = idx+n+1;
        end
        contour_pt_s(row).cross_section=row;
        contour_pt_s(row).coord=permute(contour_pt,[2 1]);
    end
    fprintf('\n')
    disp( 'All contours found. Now searching the centroids in 2D.')
    
    
    %% Recalculate centroid for 2D cross sections
    %seg_slice2D_temp=zeros(slice2D_size*2+1,slice2D_size*2+1,size(F3,1));
    c_2D=zeros(size(F3,1)-4,2);
    for row=3:size(F3,1)-2
        L=bwlabel(seg_slice2D_temp(:,:,row));
        nber_regions=max(max(L));
        emptycenter=0;
        for label=1:1:nber_regions
            [r, c] = find(L==label);
            if any(r==slice2D_size+1) && any(c==slice2D_size+1)
                break
            elseif label==nber_regions
                emptycenter=1;
            else
            end
            
        end
        
        % I now know which region to keep or whether there's been an error
        if emptycenter==1
            % Don't do anything cause there's been an error
        else % Get rid of all other regions
            X=(L==label);
            seg_slice2D_temp(:,:,row)=X;
            
        end
        
        %find the centroid of the updated image
        s=regionprops(seg_slice2D_temp(:,:,row),...
            'PixelIdxList', 'PixelList');
        %I_temp=slice2D_temp(:,:,row);
        I_temp=seg_slice2D_temp(:,:,row);
        idx = s.PixelIdxList;
        sum_region = sum(I_temp(idx)-170*ones(size(idx,1),1));
        x_region = s.PixelList(:, 1);
        y_region = s.PixelList(:, 2);
        
        c_2D(row,1) = sum(x_region .* (double(I_temp(idx))...
            -170*ones(size(idx,1),1))) / sum_region;
        c_2D(row,2) = sum(y_region .* (double(I_temp(idx))...
            -170*ones(size(idx,1),1))) / sum_region;
        
    end
    %{
    figure(9)
    for row=3:size(F3,1)-2
        subplot(round((size(F3,1))/8+1),8,row)
        imagesc(seg_slice2D_temp(:,:,row))
    end
    pause
    %}
    close all
    disp('All centroids have been found. Now calculating enclosed areas.')
    
    %% Connect edge points to form a spline perimeter
    perimeter=struct('cross_section',0,'coord',zeros(size(contour_pt)));
    %figure (13)
    
    for row=3:size(F3,1)-2
        
        perimeter(row-2).cross_section=row;
        perimeter(row-2).coord=[contour_pt_s(row).coord(:,:)];
        
        %     subplot(round((size(F3,1))/8+1),8,row)
        %     imagesc(slice2D_temp(:,:,row))
        %     hold on
        %     plot(contour_pt_s(row).coord(:,1),contour_pt_s(row).coord(:,2),...
        %         'Marker','.','LineStyle','none','MarkerEdgeColor','r','MarkerFaceColor','r')
        %     plot(c_2D(row,1),c_2D(row,2),'*')
    end
    
    
    %% Calculate the enclosed surface
    
    Areas=zeros(size(F3,1)-4,1);
    for row=3:size(F3,1)-2
        Areas(row-2)=polyarea(perimeter(row-2).coord(:,1),...
            perimeter(row-2).coord(:,2));
        
    end
    
    %% Calculate the perimeter
    P_step= 35; % number of step to interpolated around the perimeter/contour
    disp(['Now calculating perimeters using ' num2str(P_step) ' interpolated points'])
    Ps=zeros(size(F3,1)-4,1);
    perimeteri=struct('cross_section',0,'coord',zeros(P_step,2));
    fprintf([ '# of slices(-4): ' num2str(size(F3,1)-4) ...
        '. Perimeters found: '])
    for row=3:size(F3,1)-2
        if (row-3)<=10
            fprintf('\b%d',(row-3)); % update counter display
        else
            fprintf('\b\b%d',(row-3)); % update counter display
        end
        perimeteri(row-2).cross_section=row;
        perimeteri(row-2).coord=interparc(P_step, perimeter(row-2).coord(:,1),...
            perimeter(row-2).coord(:,2));
        Ps (row-2)= P_step*norm( ...
            [perimeteri(row-2).coord(3,1)-perimeteri(row-2).coord(2,1),...
            perimeteri(row-2).coord(3,2)-perimeteri(row-2).coord(2,2)]);
    end
    fprintf('\n')
    
    %% Plot interpolated perimeter
    %{
            figure (14)
            for row=3:size(F3,1)-2
                subplot(round((size(F3,1))/8+1),8,row)
                imagesc(slice2D_temp(:,:,row))
                hold on
                %plot(edge_coord(row-2).coord(:,1),edge_coord(row-2).coord(:,2),...
                %    'Marker','.','LineStyle','none','MarkerEdgeColor','r','MarkerFaceColor','r')
                plot(c_2D(row,1),c_2D(row,2),'*')
                plot(perimeteri(row-2).coord(:,1),perimeteri(row-2).coord(:,2),'r')
            end
    %}
    disp(['Now calculating the radius from ' num2str(P_step) ' interpolated points'])
    %% Recalculate the radius
    Stored_radii_para2=[];
    stored_vectors2=zeros(P_step,2);
    for row=3:1:size(F3)-2
        stored_vetors2=perimeteri(row-2).coord-...
            [c_2D(row,1)*ones(P_step,1),c_2D(row,2)*ones(P_step,1)];
        stored_radii2=zeros(size(stored_vectors2,1),1);
        for point=1:1:P_step
            stored_radii2(point)=norm(stored_vetors2(point,:));
        end
        
        Stored_radii_para2= cat(1,Stored_radii_para2,...
            [mean(stored_radii2), max(stored_radii2),min(stored_radii2)]);
    end 


end

