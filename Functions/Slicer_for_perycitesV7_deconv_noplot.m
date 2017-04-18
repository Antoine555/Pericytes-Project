% High end script that stacks the vessels images, interpolates in the ROI,
% segments with level set, skeletonises,smoothens the skeleton,
% extracts 2d slices, re calculates the centre of mass
% and then finds and plots the radius, peri, area, curvat/torsion, health,
% distance to bifurcation and exports all to excel
% this is done for all pericytes!
% if this works then I'm officially a genius...

close all
clear all


%set parameters
pause on
window=180;
slice2D_size=15;
step_size=2;
nber_slices_max= 150;

%% Set up template, import pericytes locations from excel and export their health status

% Opens the data analysis template
Template='C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Capillaries work\Functions\Results layoutV2';
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


perixyzM=xlsread(File,'sheet3','J2:L29');
%perixyzM=[199 183 29-4; 263 658 29-4];
temp=perixyzM(:,2);
perixyzM(:,2)=perixyzM(:,1);
perixyzM(:,1)=temp;
delete temp


Excel_Data=zeros(nber_slices_max-4,size(perixyzM,1)*5+1);
Excel_Data2=zeros(size(perixyzM,1),4);
Excel_Data3=zeros(size(perixyzM,1),1);

%% Stack non interpolated data
%First the deconvolved vessels
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Capillaries work\Images\vessels_deconv\26_deconv3 complete')
Stack_vessels = [];     % Set empty matrix
%t=50/256;    % Set threshold
for i=0:1:9
    Stack_vessels = cat(3,Stack_vessels, imread(sprintf('RL3 complete of C3-260%d.tif',i)));
end

for i=10:1:41
    Stack_vessels = cat(3,Stack_vessels, imread(sprintf('RL3 complete of C3-26%d.tif',i)));
end
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\3D median filter')
%Stack_vessels= ordfilt3D(Stack_vessels,14);

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

%Finally the original vessels the will be used for setting initial surface
%Second the red channel
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Capillaries work\Images')
Vref00 = [];     % Set empty matrix
%t=50/256;    % Set threshold
for i=1:1:9
    Vref00 = cat(3,Vref00, imread(sprintf('26_z00%d_c003.tif',i)));
end

for i=10:1:42
    Vref00 = cat(3,Vref00, imread(sprintf('26_z0%d_c003.tif',i)));
end

%% Start the Analysis

excluded=[];
% Start the for loop to go through one pericyte at a time
for peri=6:1:12       %size(perixyzM,1) % analalyse one pericyte at a time
    
    %store the current pericytes location to use
    perixyz=perixyzM(peri,:);
    perixyz(3)=round(perixyz(3)); %NB: this is the slice number!!!
    perixyz(1)=round(perixyz(1));
    perixyz(2)=round(perixyz(2));
    
    %% Take a volume around the pericyte
    
    %check how close to the edge the pericyte is and act upon
    if  perixyz(3)<=5 || perixyz(3)>=(42-5)  ...
            || perixyz(1)<=20 || perixyz(1)>=(size(Stack_vessels,1)-20) ...
            || perixyz(2)<=20 || perixyz(2)>=(size(Stack_vessels,2)-20)
        % don't analyse this pericyte and go to the next
        excluded=cat(1,excluded,peri); % record which have been excluded
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
    
    %% interpolate this volume and call it V
    
    V0_db=double(V0);
    ny=size(V0,2);nx=size(V0,1);nz=size(V0,3)*4.817734273; %% desired output dimensions
    [yq,xq,zq]=  ndgrid(linspace(1,size(V0,2),ny),...
        linspace(1,size(V0,1),nx),...
        linspace(1,size(V0,3),nz));
    V=interp3(V0_db,yq,xq,zq);
    
    
    Vref0_db=double(Vref0);
    Vref=interp3(Vref0_db,yq,xq,zq);
    
    
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
        
        pause
    %}
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
        
        pause
    %}
    
    %% Cleaning and labelling the skeleton
    % Now use Skel2Graph on the skeleton
    addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Skel2graph3D')
    
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
        
        pause
    %}
    %% Continue the analyisis if the link found is correct
    %skip this pericyte if the correct link could not be computed
    if minDistance>=20 % link is too far from perictyte
        excluded=cat(1,excluded,peri);
        disp(['Pericyte' num2str(peri) 'excluded'])
    else
        
        disp(['measuring initiated for pericyte' num2str(peri)])
        
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
        %if minDistance < overallMinDistance
        % This boundary is the closer one.
        %end
        %closest_ptcoord=F3(indexOfMin,:);
        
        
        %% Extract Slices
        
        %figure (7)
        % store value for the latest skeleton
        c_skel4=[x_skel4,y_skel4,z_skel4];
        %pause on
        %map=[0 0 0; 0 0.2 0;0 0.4 0;0 0.6 0; 0 0.8 0;0 1 0];
        
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
            
            % eaxtract slices if using the point on smoot spline as origin:
            [slice2D_temp(:,:,row),sliceInd_temp(:,:,row),...
                subX_temp(:,:,row),subY_temp(:,:,row),subZ_temp(:,:,row)]...
                =extractSlice(V,...
                F3(row,1),F3(row,2),F3(row,3),d(1),d(2),d(3),slice2D_size);
            
            
            
        end
        slice2D_temp(isnan(slice2D_temp)) = 0 ;
        %{
            for row=1:size(F3,1)
                %Plot slices as we circulate through the centrepoints
                surf(subX_temp(:,:,row),subY_temp(:,:,row),subZ_temp(:,:,row),...
                    slice2D_temp(:,:,row),'FaceColor','texturemap','EdgeColor','none');
                axis([1 S2(1) 1 S2(2) 1 S2(3)]); colormap(map);
                view(-17,46);
                drawnow
                hold on
            end
            
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
            
            % Plot all the 2D slices
            figure(8)
            for row=1:size(F3,1)
                subplot(round((size(F3,1))/8+1),8,row)
                imagesc(slice2D_temp(:,:,row))
            end
            
            pause
        %}
        
        %% find the contour using hybrid level set
        
        addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Level set toolbox\AOSLevelsetSegmentationToolboxM')
        seg_slice2D_temp=zeros(slice2D_size*2+1,slice2D_size*2+1,size(F3,1));
        contour_ptM=[];
        phi=[];
        contour_pt_s=struct('cross_section',0,'coord',[]);
        for row=1:1:size(F3,1)
            slice2D_temp(:,:,row)=slice2D_temp(:,:,row)./max(max(slice2D_temp(:,:,row)))*255;
            %I =  ordfilt3D(slice2D_temp(:,:,row),14);
            I =  slice2D_temp(:,:,row);
            %phi = ac_SDF_2D('circle', size(I), c_2D(15,:),slice2D_size) ;
            phi = ac_SDF_2D('circle', size(I), [round(size(I,1)/2),round(size(I)/2)] ,2);
            
            g = ac_gradient_map(I, 5);
            mu = 0.1*slice2D_temp(16,16,row);
            propagation_weight = .6; GAC_weight = 5000;
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
        
        
        %% Recalculate centroid for 2D cross sections
        %seg_slice2D_temp=zeros(slice2D_size*2+1,slice2D_size*2+1,size(F3,1));
        c_2D=zeros(size(F3,1),2);
        for row=1:size(F3,1)
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
            for row=1:size(F3,1)
                subplot(round((size(F3,1))/8+1),8,row)
                imagesc(seg_slice2D_temp(:,:,row))
            end
            
            pause
            
        %}
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
        
        % pause
        
        %% Calculate the enclosed surface
        
        Areas=zeros(size(F3,1)-4,1);
        for row=3:size(F3,1)-2
            Areas(row-2)=polyarea(perimeter(row-2).coord(:,1),...
                perimeter(row-2).coord(:,2));
            
        end
        
        %% Calculate the perimeter
        Ps=zeros(size(F3,1)-4,1);
        perimeteri=struct('cross_section',0,'coord',zeros(40,2));
        for row=3:size(F3,1)-2
            perimeteri(row-2).cross_section=row;
            perimeteri(row-2).coord=interparc(40, perimeter(row-2).coord(:,1),...
                perimeter(row-2).coord(:,2));
            Ps (row-2)= 40*norm( ...
                [perimeteri(row-2).coord(3,1)-perimeteri(row-2).coord(2,1),...
                perimeteri(row-2).coord(3,2)-perimeteri(row-2).coord(2,2)]);
        end
        
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
            
            pause
        %}
        %% Recalculate the radius
        Stored_radii_para2=[];
        stored_vectors2=zeros(40,2);
        for row=3:1:size(F3)-2
            stored_vetors2=perimeteri(row-2).coord-...
                [c_2D(row,1)*ones(40,1),c_2D(row,2)*ones(40,1)];
            stored_radii2=zeros(size(stored_vectors2,1),1);
            for point=1:1:40
                stored_radii2(point)=norm(stored_vetors2(point,:));
            end
            
            Stored_radii_para2= cat(1,Stored_radii_para2,...
                [mean(stored_radii2), max(stored_radii2),min(stored_radii2)]);
        end
        
        %% Plot radius, surface and perimeter along arclength
        figure %(15)
        row=3:1:size(F3,1)-2;
        
        subplot(2,2,1)
        %radius
        plot(row, Stored_radii_para2(:,1),'r')    %plot mean
        hold on
        plot(row, Stored_radii_para2(:,2),'g')    %plot max
        plot(row, Stored_radii_para2(:,3))        %plot min in blue
        hold on
        plot (indexOfMin, Stored_radii_para2(indexOfMin-2,1), 'Marker','.','MarkerSize',8,...
            'MarkerFaceColor','g','MarkerEdgeColor','m')
        
        %Area
        subplot(2,2,2)
        plot(row,Areas,'m')
        hold on
        plot (indexOfMin, Areas(indexOfMin-2), 'Marker','.','MarkerSize',8,...
            'MarkerFaceColor','g','MarkerEdgeColor','g')
        
        %Perimeter
        subplot(2,2,3)
        plot(row,Ps,'c')
        hold on
        plot (indexOfMin, Ps(indexOfMin-2), 'Marker','.','MarkerSize',8,...
            'MarkerFaceColor','g','MarkerEdgeColor','m')
        
        %Area/perimeter and max= rmean/2 (for a circle)
        shape_para=Areas./Ps;
        subplot(2,2,4)
        plot(row,shape_para,'g',row, Stored_radii_para2(:,1)/2,'r' )
        hold on
        plot (indexOfMin, shape_para(indexOfMin-2), 'Marker','.','MarkerSize',8,...
            'MarkerFaceColor','g','MarkerEdgeColor','m')
        
        %pause
        
        %% make a 2D projection
        
        %fig=figure (16);
        figure
        MIP= imread('26Composite.tif (RGB).tif');
        imshow(MIP(xe1:xe2,...
            ye1:ye2,:));
        hold on
        plot (mylink_x,mylink_y,'.')
        plot(x_peri_wind,y_peri_wind,'o', 'Color', 'y','MarkerFaceColor','y')
        plot(F3(indexOfMin,1),F3(indexOfMin,2),'Marker','+','MarkerSize',10,...
            'MarkerFaceColor','r','MarkerEdgeColor','r')
        
        drawnow
        %print(fig,'-dpng');
        
        %pause
        
        %% Plot for report
        %{
            figure (17)
            row=3:1:size(F3,1)-2;
            
            %radius
            plot(row(1:end), Stored_radii_para2((1:end),1)*0.415,'r')    %plot mean
            hold on
            plot(row(1:end), Stored_radii_para2((1:end),2)*0.415,'g')    %plot max
            plot(row(1:end), Stored_radii_para2((1:end),3)*0.415)        %plot min in blue
            %hold on
            %plot (indexOfMin, Stored_radii_para2(indexOfMin-2,1), 'Marker','.','MarkerSize',8,...
            %    'MarkerFaceColor','g','MarkerEdgeColor','m')
            figure
            %Area
            [AX,H1,H2]=plotyy(row(1:end),Areas(1:end)*(0.4151^2),row(1:end),Ps(1:end)/2.5*0.4151);
            %hold on
            %plot (indexOfMin, Areas(indexOfMin-2), 'Marker','.','MarkerSize',8,...
            %     'MarkerFaceColor','g','MarkerEdgeColor','g')
            xlabel('Slice number')
            
            ylabel(AX(1),'Area (micrometers^2)') % left y-axis
            ylabel(AX(2),'Perimeter (micrometers)') % right y-axis
            
            pause
        %}
        %% Frenet-Serret
        %{
            addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Frenet')
            %[F4] = interparc(nber_slices*5,F(:,1),F(:,2),F(:,3));
            kappalist=[];
            taulist=[];
            ttlist=[];
            nnlist=[];
            bblist=[];
            score=[];
            [kappalist,taulist,ttlist,nnlist,bblist,score]=frenet_robust(F2.',25,0.2);
            %if nanmean(kappalist0.*kappalist')<0
            %    kappalist=-kappalist;
            %end;
            indFS=1:size(F2);
            indFS=indFS/5;
            % plot
            figure
            [AX,H1,H2]=plotyy(indFS(1:150),kappalist(1:150),indFS(1:150),taulist(1:150));
            %hold on
            %plot (indexOfMin, Areas(indexOfMin-2), 'Marker','.','MarkerSize',8,...
            %     'MarkerFaceColor','g','MarkerEdgeColor','g')
            xlabel('Slice number')
            
            ylabel(AX(1),'Curvature') % left y-axis
            ylabel(AX(2),'Torsion') % right y-axi
        %}
        
        %% plot pericyte locations
        
        %figure;
        %MIP= imread('26Composite.tif (RGB).tif');
        %imshow(MIP(435-window/2:435+window/2,...
        %    422-window/2:422+window/2,:));
        %hold on
        
        %plot(window/2+1,window/2+1,'o', 'Color', 'y','MarkerFaceColor','y')
        
        Excel_Data(1:(size(F3,1)-4),1)=(3:1:size(F3,1)-2);
        Excel_Data(1:(size(F3,1)-4),(peri*5-3):(peri*5-1))=Stored_radii_para2;
        Excel_Data(1:(size(F3,1)-4),(peri*5))=Areas;
        Excel_Data(1:(size(F3,1)-4),(peri*5+1))=Ps;
        %Excel_Data(1:(size(F3,1)-4),(peri*5+2))=kappalist(1:26);
        %Excel_Data(1:(size(F3,1)-4),(peri*5+3))=taulist(1:26);
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
end


xlswrite(File,Excel_Data,'sheet2','A3');
xlswrite(File,Excel_Data3,'sheet1','E3');
xlswrite(File,Excel_Data2,'sheet1','F3');