% High end script that stacks the vessels images, interpolates in the ROI,
% segments with level set, skeletonises, smoothens the skeleton,
% extracts 2d slices, re calculates the centre of mass
% and then finds and plots the radius, peri, area, curvat/torsion, health,
% distance to bifurcation and exports all to excel
% this is done for all pericytes!
% for this the window is defined differently when the pericyte is close
% to the edge of the volume.
% V9 uses the updated (Jan 2016) version of Skel2Graph3D which keeps the
% endpoints in the skeleton
% xy version will simply measure along a z-plane at the centroid

%close all
%clear all


%set parameters
pause on
window=180;
slice2D_size=15;
step_size=2;
nber_slices_max= 150;
aspect_ratio=4.817734273;

%% Set up template, import pericytes locations from excel and export their health status

% Opens the data analysis template
Template='C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Capillaries work\Functions\Results layout xy';
% Initiates Excel COM Object
Excel = actxserver ('Excel.Application');
set(Excel, 'Visible', 1);
template=invoke(Excel.Workbooks,'Open',Template);
% Allows user to select output file
Exportdir=uigetdir('','Please select export directory');
prompt = {'Please enter EXCEL File name:'};
name = 'Data Output';
numlines = 1;
acq_date=date;
defaultanswer = {['stack nb 26','_',acq_date,'.xlsx']};
options.Resize = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
FileName = cell2mat(inputdlg(prompt,name,numlines,defaultanswer,options));
File=fullfile(Exportdir,FileName);
% Saves the template to the user selected output file
invoke(template,'SaveCopyAs',File);
invoke(template, 'Close');
% Quit Excel
invoke(Excel, 'Quit');
% End process
delete(Excel);


perixyzM=xlsread(File,'sheet3','J2:L50');
%perixyzM=[199 183 29-4; 263 658 29-4];
temp=perixyzM(:,2);
perixyzM(:,2)=perixyzM(:,1);
perixyzM(:,1)=temp;
%delete temp


Excel_Data=zeros(nber_slices_max-4,size(perixyzM,1)*5+1);
Excel_Data2=zeros(size(perixyzM,1),4);
Excel_Data3=zeros(size(perixyzM,1),1);

%% Stack non interpolated data
%First the deconvolved vessels
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Capillaries work\Images\vessels_deconv\26_deconv3 complete')
Stack_vessels = [];     % Set empty matrix
%t=50/256;    % Set threshold
disp('Stacking started')
for i=0:1:9
    Stack_vessels = cat(3,Stack_vessels, imread(sprintf('RL3 complete of C3-260%d.tif',i)));
end

for i=10:1:41
    Stack_vessels = cat(3,Stack_vessels, imread(sprintf('RL3 complete of C3-26%d.tif',i)));
end
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\3D median filter')
%Stack_vessels= ordfilt3D(Stack_vessels,14);
disp('deconvolved green/vessel channel stacked')
%Second the red channel
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Capillaries work\Images')
V_red = [];     % Set empty matrix
%t=50/256;    % Set threshold
for i=1:1:9
    V_red = cat(3,V_red, imread(sprintf('26_z00%d_c002.tif',i)));
end

for i=10:1:42
    V_red = cat(3,V_red, imread(sprintf('26_z0%d_c002.tif',i)));
end
disp('red channel stacked')

%Finally the original vessels the will be used for setting initial surface
% and in this script to find the maximum radius
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Capillaries work\Images')
Vref00 = [];     % Set empty matrix
%t=50/256;    % Set threshold
for i=1:1:9
    Vref00 = cat(3,Vref00, imread(sprintf('26_z00%d_c003.tif',i)));
end

for i=10:1:42
    Vref00 = cat(3,Vref00, imread(sprintf('26_z0%d_c003.tif',i)));
end
disp('undeconvolved green/vessel channel stacked to intialise level set')

%% Find the look up table for the interpolated z-coordinates

[ InterpTable ] = InterpLookUp( size(Vref00,3), aspect_ratio );

%% Start the Analysis

excluded=[];
% Start the for loop to go through one pericyte at a time
for peri=1:1:size(perixyzM,1)       %size(perixyzM,1) % analalyse one pericyte at a time
    disp(['pericyte ' num2str(peri) ' selected.'])
    %store the current pericytes location to use
    perixyz=perixyzM(peri,:);
    perixyz(3)=round(perixyz(3)); %NB: this is the slice number!!!
    perixyz(1)=round(perixyz(1));
    perixyz(2)=round(perixyz(2));
    
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
    %zi=(perixyz(3)-1)*4.818-1;  % This is the coordinate in the interpolated volume!!
    % Look up the interpolated z-coordinate
    if ~isinteger(perixyz(3));
        zi=(InterpTable(floor(perixyz(3)),2)+InterpTable(ceil(perixyz(3)),2))/2;
    else
    zi=InterpTable(perixyz(3),2);
    end
    
    disp( ['Pericyte is valid and volume window has been selected.'...
        ' Interpolation of this volume started.'])
    
    %% interpolate this volume (deconv and raw) and call it V
    
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
    
    
    % initialise parameters
    smooth_weight = 2;
    image_weight = 1e-4;
    delta_t = 1;
    margin = 10;
    %phi = zeros(size(V)); % alternate initialisation
    %phi(margin:end-margin, margin:end-margin, margin:end-margin) = 1;
    phi=Vref>130;
    phi = ac_reinit(phi-.5);
    
    % run ChanVese level set on deconvolved volume
    for i = 1:2
        phi = ac_ChanVese_model(V, phi, smooth_weight, image_weight, delta_t, 1);
        
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
        subplot(2,4,i); imshow(squeeze(V(:,slice(i),:)),[]); hold on;
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
    imOut1=V;
    
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
    subplot(2,2,1); imagesc(Vref(:,:,25))
    subplot(2,2,2); imagesc(squeeze(Vref(:,50,:)))
    subplot(2,2,3); imagesc(phi_seg(:,:,25))
    subplot(2,2,4); imagesc(squeeze(phi_seg(50,:,:)))
    drawnow
    
    
    % re-initialise parameters for level set of Vref
    smooth_weight = 2;
    image_weight = 1e-4;
    delta_t = 1;
    margin = 10;
    %phi = zeros(size(V)); % alternate initialisation
    %phi(margin:end-margin, margin:end-margin, margin:end-margin) = 1;
    phi2=Vref>130;
    phi2 = ac_reinit(phi2-.5);
    
    % run ChanVese level set on deconvolved volume
    for i = 1:2
        phi2 = ac_ChanVese_model(Vref, phi2, smooth_weight, image_weight, delta_t, 1);
        
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
        subplot(2,4,i); imshow(squeeze(Vref(:,slice(i),:)),[]); hold on;
        c2 = contours(squeeze(phi2(:,slice(i),:)),[0,0]);
        zy_plot_contours(c2,'linewidth',2);
    end
    %find contour of deconvolved volume
    phi_seg2=phi2>0;
    
    figure;
    slice = [30,35,40,45,50,55,60,65];
    for i = 1:8
        subplot(2,4,i); imshow(squeeze(phi_seg2(:,slice(i),:)),[]); hold on;
        c2 = contours(squeeze(phi2(:,slice(i),:)),[0,0]);
        zy_plot_contours(c2,'linewidth',2);
    end
    
    
    seg_imOut2=(phi_seg2);
    
    
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
    subplot(2,2,1); imagesc(Vref(:,:,25))
    subplot(2,2,2); imagesc(squeeze(Vref(:,50,:)))
    subplot(2,2,3); imagesc(phi_seg2(:,:,25))
    subplot(2,2,4); imagesc(squeeze(phi_seg2(50,:,:)))
    drawnow
    disp('Segmentation finished. Skeletonisation starting')
    
    
    %% Skeletonisation
    
    addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Skeleton3D')
    skel = Skeleton3D(seg_imOut);
    
    figure
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
    
    w=size(skel,1);
    l=size(skel,2);
    h=size(skel,3);
    [x_skel,y_skel,z_skel]=ind2sub([w,l,h],find(skel(:)));
    plot3(x_skel,y_skel,z_skel,'.','Markersize',8,'MarkerFaceColor','b','Color','b');
    set(gcf,'Color','White');
    view(3)
    drawnow
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
    
    % display result
    figure;
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
    figure
    w=size(skel3,1);
    l=size(skel3,2);
    h=size(skel3,3);
    skel5=skel3-skel4;
    [x_skel5,y_skel5,z_skel5]=ind2sub([w,l,h],find(skel5(:)));
    plot3(x_skel5,y_skel5,z_skel5,'.','Markersize',8,'MarkerFaceColor','b','Color','b');
    axis ([0 S2(1) 0 S2(2) 0 S2(3)]);
    drawnow
    hold on
    
    w=size(skel4,1);
    l=size(skel4,2);
    h=size(skel4,3);
    [x_skel4,y_skel4,z_skel4]=ind2sub([w,l,h],find(skel4(:)));
    plot3(x_skel4,y_skel4,z_skel4,'.','Markersize',8,'MarkerFaceColor','r','Color','r');
    view(3)
    drawnow
    hold on
    axis tight; axis off;
    daspect([1,1,1])
    view(3);
    
    plot3 (x_peri_wind,y_peri_wind,zi,'MarkerSize',7,'MarkerFaceColor','m','MarkerEdgeColor','m','Marker','o')
    
    
    
    %% Continue the analyisis if the link found is correct
    %skip this pericyte if the correct link could not be computed
    if minDistance>=20 % link is too far from perictyte
        excluded=cat(1,excluded,peri);
        disp(['Pericyte' num2str(peri) 'excluded as too far from'...
            ' nearest link. Check z coordinate or increase window size.'...
            ' Closing open figures then going go to next peri.'])
        close all
        pause (2)
        continue
    else
    end
    
    
    %% Use spline interpolation to smooth the centreline
    
    addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\SmoothN');
    figure 
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
    plot3 (mylink_x,mylink_y,mylink_z,'.',u{1},u{2},u{3},'r')
    axis tight; axis off;
    daspect([1,1,1])
    view(3)
    box on
    disp('Link found and smoothed. Check plots then press any key.')
    pause (1)
    
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
            ' analyse these. Closing all open windows then' ...
            ' and going to the next pericyte.'])
        close all
        pause (2)
        continue
    else
    end
    disp('Starting taking measurements along x-y plane only.')
    
    %% Extract Slices
    %% find the contour using hybrid level set
    %% Recalculate centroid for 2D cross sections
    %% Connect edge points to form a spline perimeter
    %% Calculate the enclosed surface
    %% Calculate the perimeter
    %% Plot interpolated perimeter
    %% Recalculate the radius
    
    %% Measure the diameter lying on the x-y plane
    xyradius1=zeros((size(F3,1)-4)); % matrix for recording the two point
    xyradius2=zeros((size(F3,1)-4)); % matrix for recording the two point
    v= [1,0,0];%initialise direction for ray casting
    diameter=zeros((size(F3,1)-4),1);
    % Find the radii for each centrepopint
    addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Timed Progress Bar')
    targetWorkCount = size(F3,1)-2;
    barWidth= int32( targetWorkCount/3 );
    p =  TimedProgressBar( targetWorkCount, barWidth, ...
    ' measuring, wait for ', ', completed ', 'Concluded in ' );
                         % Also percent = p.stop;
    parfor row=3:1:size(F3,1)-2
        
        % Find direction vector using derivative at point
        d=dF3dt(row,:);
        d=d/max(d)*50;
        if d(1)==0 && d(2)==0 % cross section is on a z-plane
            v= [1,0,0];
            % take the direction of previous step v=v
        else
            v=[-d(2), d(1), 0]/(sqrt(d(2)^2+d(1)^2));  % find v
        end
        %first direction
        vr= F3(row,:); %Start at the centroid
        vr_value=interp3(phi2, vr(2),vr(1),vr(3));
        e=0.1; %set tolerance on accuracy of radius
        vr_value_old=vr_value;
        while vr_value >=0 %cast ray and check value
            vr_value_old=vr_value; % save the before last value
            vr=vr+v;
            vr_value= interp3(phi2, vr(2),vr(1),vr(3));
        end
        % Get as close as possible to the point where vr_value=0
        vr3=(vr*abs(vr_value)+(vr-v)*vr_value_old)/(abs(vr_value)+vr_value_old);
        vr3_value=interp3(phi2, vr3(2),vr3(1),vr3(3));
        while abs(vr3_value)>e
            vr3=vr3+0.1*v*vr_value/abs(vr_value);
            vr3_value=interp3(phi2, vr3(2),vr3(1),vr3(3));
        end
        xyradius1(row-2)=norm(vr3-F3(row,:)); % record the corresponding rad
        
        % opposite direction
        vr= F3(row,:); %Start at the centroid
        vr_value=interp3(phi2, vr(2),vr(1),vr(3));
        v=-v;
        vr_value_old=vr_value;
        while vr_value >=0 %cast ray and check value
            vr_value_old=vr_value; % save the before last value
            vr=vr+v;
            vr_value= interp3(phi2, vr(2),vr(1),vr(3));
        end
        % Get as close as possible to the point where vr_value=0
        vr3=(vr*abs(vr_value)+(vr-v)*vr_value_old)/(abs(vr_value)+vr_value_old);
        vr3_value=interp3(phi2, vr3(2),vr3(1),vr3(3));
        while abs(vr3_value)>e
            vr3=vr3+0.1*v*vr_value/abs(vr_value);
            vr3_value=interp3(phi2, vr3(2),vr3(1),vr3(3));
        end
        xyradius2(row-2)=norm(vr3-F3(row,:)); % record the corresponding rad
        diameter(row-2)= xyradius1(row-2)+xyradius2(row-2);
        p.progress;
    end
    p.stop;  
    disp( 'Finished measuring.')
    
    %% PloT diameter
    figure
    row=3:1:size(F3,1)-2;
    %radius
    plot(row(1:end), diameter*0.415,'g')    %plot mean
    hold on
    plot (indexOfMin, diameter(indexOfMin-2)*0.415, 'Marker','.','MarkerSize',8,...
        'MarkerFaceColor','g','MarkerEdgeColor','m')
    
    
    %% make a 2D projection
    
    fig=figure (16);
    MIP= imread('26Composite.tif (RGB).tif');
    imshow(MIP(xe1:xe2,...
        ye1:ye2,:));
    hold on
    plot (mylink_x,mylink_y,'.')
    plot(x_peri_wind,y_peri_wind,'o', 'Color', 'y','MarkerFaceColor','y')
    plot(F3(indexOfMin,1),F3(indexOfMin,2),'Marker','+','MarkerSize',10,...
        'MarkerFaceColor','r','MarkerEdgeColor','r')
    
    print(fig,'-dpng');
    disp('Closing all open figures and going to next pericyte.')
    pause (3)
    close all 
    pause (1)
    
    %% Plot for report
    %% Frenet-Serret
    %% plot pericyte locations
    
    %% Save measurements data inside a matrix
    Excel_Data(1:(size(F3,1)-4),1)=(3:1:size(F3,1)-2);
    Excel_Data(1:(size(F3,1)-4), peri+1)=diameter;
    Excel_Data2(peri,:)=[indexOfMin, indexOfMin*0.4151*2, mylink_al, mylink_al/mylink_l];
    
    healthM= V_red(round(perixyz(1))-3 : round(perixyz(1))+3,...
        round(perixyz(2))-3: round(perixyz(2)+3),...
        (round(perixyz(3)/4.812+1)-1):(round(perixyz(3)/4.812+1)));
    health =sum(sum(sum(healthM)))/(7*7*3);
    if health>=20
        Excel_Data3(peri)=1;
    else
        Excel_Data3(peri)=0;
    end
end

disp('All pericyte have been explored and data is being exported to Excel.')

xlswrite(File,Excel_Data,'sheet2','A3');
xlswrite(File,Excel_Data3,'sheet1','E3');
xlswrite(File,Excel_Data2,'sheet1','F3');

disp('Script finished!')