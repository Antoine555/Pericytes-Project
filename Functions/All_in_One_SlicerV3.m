% High end script that stacks the vessels images, segment with graph cut or
% level set method, skeletonises, soothen the skeleton, extract 2d slices,
% calculate the centre of mass and then finds ond plaot the radius evolution.

window=150;
slice2D_size=15;
nber_slices= 60;
peryxyz=[110 110 29-4];
   V=V0(peryxyz(1)-window/2:peryxyz(1)+window/2,peryxyz(2)-window/2:peryxyz(2)+window/2,:);
    Vref=Vref0(peryxyz(1)-window/2:peryxyz(1)+window/2,peryxyz(2)-window/2:peryxyz(2)+window/2,:);
    S=size(V);
        x=linspace (1,S(1),S(1));
    y=linspace(1,S(2),S(2));
    z=linspace (1,S(3),S(3));

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
        
        if exist('h','var') && all(ishandle(h)), delete(h); end
        iso = isosurface(phi);
        h = patch(iso,'facecolor','w');  axis equal;  view(3);
        set(gcf,'name', sprintf('#iters = %d',i));
        view(3);
        drawnow;
    end
    
    %%
    figure;
    slice = [30,35,40,45,50,55,60,65];
    for i = 1:8
        subplot(2,4,i); imshow(squeeze(V(:,slice(i),:)),[]); hold on;
        c = contours(squeeze(phi(:,slice(i),:)),[0,0]);
        zy_plot_contours(c,'linewidth',2);
    end
    
    phi_seg=phi>0;
    
    figure;
    slice = [30,35,40,45,50,55,60,65];
    for i = 1:8
        subplot(2,4,i); imshow(squeeze(phi_seg(:,slice(i),:)),[]); hold on;
        c = contours(squeeze(phi(:,slice(i),:)),[0,0]);
        zy_plot_contours(c,'linewidth',2);
    end
    
    
    seg_imOut=(phi_seg);
    imOut1=(Vref);
    
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
    
    keyboard
      %check that you are happy with the segmentation

%% Skeletonisation

addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Skeleton3D')
skel = Skeleton3D(permute(seg_imOut,[2 1 3]));

%addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Fast Marching')
%skel=skeleton(seg_stack3);

% Plotting
figure(3)
col=[.7 .7 .8];
hiso = patch(isosurface(seg_imOut,0),'FaceColor',col,'EdgeColor','none');
hiso2 = patch(isocaps(seg_imOut,0),'FaceColor',col,'EdgeColor','none');
axis equal;axis off;
lighting phong;
isonormals(seg_imOut,hiso);
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

%% Cleaning and labelling the skeleton
% Now use Skel2Graph on the skeleton
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Skel2graph3D')

skel2=skel;

w = size(skel2,1);
l = size(skel2,2);
h = size(skel2,3);

% convert skeleton to graph structure
[A,node,link] = Skel2Graph3D(skel2,20);

% convert graph structure back to (cleaned) skeleton
skel3 = Graph2Skel3D(node,link,w,l,h);

i=1;
% iteratively convert until there are no more 2-nodes left
[A2,node2,link2] = Skel2Graph3D(skel3,20);
while(min(cellfun('length',{node2.conn}))<3)
    skel3 = Graph2Skel3D(node2,link2,w,l,h);
    [A2,node2,link2] = Skel2Graph3D(skel3,20);
    i=i+1;
    if i>20
        break
    end
end;

% display result
figure(4);
hold on;
xn=[];
yn=[];
zn=[];
labels=[];
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
hiso = patch(isosurface(seg_imOut,0),'FaceColor',col,'EdgeColor','none');
hiso2 = patch(isocaps(seg_imOut,0),'FaceColor',col,'EdgeColor','none');
axis equal;%axis off;
lighting phong;
isonormals(seg_imOut,hiso);
alpha(0.5);
%set(gca,'DataAspectRatio',[1 1 1])
camlight;

%% Extract a single link of the skeleton

%Plot the cleaned skeleton with annotated nodes
figure(5)
w=size(skel3,1);
l=size(skel3,2);
h=size(skel3,3);
[x,y,z]=ind2sub([w,l,h],find(skel3(:)));
plot3(x,y,z,'.','Markersize',8,'MarkerFaceColor','b','Color','b');
drawnow
hold on

annot_gap=6*ones(size(yn));
plot3(xn,yn,zn,'o','Markersize',9,...
    'MarkerFaceColor','y',...
    'Color','k');
text(xn+annot_gap,yn+annot_gap,zn+annot_gap,labels)
view(-17,46);

%ask the user to input two node numbers connected by the link of interest
prompt = ...
    'Which link? [smaller node number,larger node number]=';
mynodes=input(prompt);
for nodes_index=1:size(link2,2)
    
    if link2(nodes_index).n1==mynodes(1) &&  link2(nodes_index).n2==mynodes(2)
        %|| link2(nodes_index).n1=mynodes(2) && link2(nodes_index).n2=mynodes(1)
        
        mylink=link2(nodes_index).point;
        break
    end
end

% Create the new matrix holding the single link
skel4=zeros(size(skel3));
skel4(mylink)=1;

figure(6)
w=size(skel4,1);
l=size(skel4,2);
h=size(skel4,3);
[x_skel4,y_skel4,z_skel4]=ind2sub([w,l,h],find(skel4(:)));
plot3(x_skel4,y_skel4,z_skel4,'.','Markersize',8,'MarkerFaceColor','b','Color','b');

drawnow
hold on

plot3(xn(mynodes),yn(mynodes),zn(mynodes),'o','Markersize',9,...
    'MarkerFaceColor','y',...
    'Color','k');
axis ([0 S2(1) 0 S2(2) 0 S2(3)]);
view(-17,46);
drawnow

pause

%% Create a reduced 'box' volume around the link to which we apply  (will do this later)
% the Frangi filter. This help reducing computation time.

%% Frangi filter for extracting the eigenvalues

%{
addpath ('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Frangi2')
% mex eig3volume.c

% Frangi Filter the stent volume
options.BlackWhite=false;
options.FrangiScaleRange=[1 10];
options.FrangiScaleRatio=2;
options.FrangiC=10;
[imOut_filtered,Scale,Vx,Vy,Vz]=FrangiFilter3D(imOut,options);

%   % Show maximum intensity plots of input and result
figure(7),
subplot(2,2,1), imshow(squeeze(max(imOut,[],2)),[])
subplot(2,2,2), imshow(squeeze(max(imOut_filtered,[],2)),[])
subplot(2,2,3), imshow(imOut(:,:,50),[])
subplot(2,2,4), imshow(imOut_filtered(:,:,50),[])
%}

%% Use spline interpolation to smooth the centreline

figure(7)

%plot centrepoint
w=size(skel4,1);
l=size(skel4,2);
h=size(skel4,3);
[x_skel4,y_skel4,z_skel4]=ind2sub([w,l,h],find(skel4(:)));
plot3(x_skel4,y_skel4,z_skel4,'.','Markersize',8,'MarkerFaceColor','b','Color','b');
drawnow
hold on

% get coordinate of the link points only
[mylink_x,mylink_y,mylink_z]=ind2sub(size(skel4),mylink);
mylink_coord=[mylink_x;mylink_y;mylink_z];
%fnplt(cscvn(mylink_coord(:,1:end )),'r',2)

[pp]=csaps(mylink_z,[mylink_x;mylink_y]);
val=fnval(pp,mylink_z);
%plot3(mylink_x,mylink_y,mylink_z);
plot3(val(1,:),val(2,:),mylink_z,'r-')
grid on

drawnow;
view(-17,46);
hold on
col=[.7 .7 .8];
hiso = patch(isosurface(seg_imOut,0),'FaceColor',col,'EdgeColor','none');
hiso2 = patch(isocaps(seg_imOut,0),'FaceColor',col,'EdgeColor','none');
axis equal;%axis off;
lighting phong;
isonormals(seg_imOut,hiso);
alpha(0.5);
%set(gca,'DataAspectRatio',[1 1 1])
camlight;

pause

%% interpolate for equidistant points along the smooth centreline
%http://fr.mathworks.com/matlabcentral/fileexchange/34874-interparc
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\interparc')

nber_slices= 15; % number of points to be interpolated i.e.number of slices
F=permute([val;mylink_z],[2 1]); %Matrix of coordinates for the smooth spline
F2=unique(F,'rows'); %remove repeating points
[F3,dF3dt] = interparc(nber_slices,F2(:,1),F2(:,2),F2(:,3)); % also returns derivatives, and a function handle


%% Find the derivative of the piece wise polynomial

%% Radius measuring from smooth spline + 2D centre of mass

figure (8)
% store value for the latest skeleton
c_skel4=[x_skel4,y_skel4,z_skel4];
pause on
map=[0 0 0; 0 0.2 0;0 0.4 0;0 0.6 0; 0 0.8 0;0 1 0];

%Extract slices and plot them onto isosurface
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Extract Slice')
Stored_radii_para=[];%set matrix for storing mean max min radii at each centrepoint
d=[0 0 0]; %initialise direction vector
slice2D_size= 15; %choose parametric size of the slices taken

slice2D_temp=zeros(31,31,size(F3,1));
sliceInd_temp=zeros(31,31,size(F3,1));
subX_temp=zeros(31,31,size(F3,1));
subY_temp=zeros(31,31,size(F3,1));
subZ_temp=zeros(31,31,size(F3,1));
%c_2D=zeros(1,size(F3,1)-1); %centre of image is fixed for now

for row=1:size(F3,1) %for all centrepoints
    % Find direction vector (more than one method)
    % d=F3(row,:)-F3(row-1,:); % note cannot start at row one
    % d=d*5;
    
    % Find direction vector using derivative at point
    d=dF3dt(row,:);
    d=d*5;
    
    % eaxtract slices if using the point on smoot spline as origin:
    [slice2D_temp(:,:,row),sliceInd_temp(:,:,row),...
        subX_temp(:,:,row),subY_temp(:,:,row),subZ_temp(:,:,row)]...
        =extractSlice(permute(imOut1, [2 1 3]),...
        F3(row,1),F3(row,2),F3(row,3),d(1),d(2),d(3),slice2D_size);
    
end

for row=1:size(F3,1)
    %Plot slices as we circulate through the centrepoints
    surf(subX_temp(:,:,row),subY_temp(:,:,row),subZ_temp(:,:,row),slice2D_temp(:,:,row),'FaceColor','texturemap','EdgeColor','none');
    axis([1 S2(1) 1 S2(2) 1 S2(3)]); colormap(map);
    view(-17,46);
    drawnow
    hold on
end

pause
col=[.7 .7 .8];
hiso = patch(isosurface(seg_imOut,0),'FaceColor',col,'EdgeColor','none');
hiso2 = patch(isocaps(seg_imOut,0),'FaceColor',col,'EdgeColor','none');
axis equal;%axis off;
lighting phong;
isonormals(seg_imOut,hiso);
alpha(0.5);
%set(gca,'DataAspectRatio',[1 1 1])
camlight;
pause

% Plot all the 2D slices
figure(9)
for row=1:size(F3,1) 
subplot(round((size(F3,1))/8+1),8,row)
imagesc(slice2D_temp(:,:,row))
end

pause

% Measure the radii from these images
for row=1:1:size(F3,1) %for all centrepoints
    nber_radii=4; %choose the number of radii to be measured per centrepoint
    angle_incr=360/nber_radii;
    stored_radii=zeros(1,nber_radii); %set an empty matrix for storing measured radii
    c_2D=slice2D_size*2/2+[1 1];   %set coordinates of the centrepoint
    flag=0;  %set a flag to be able to break out of a nested while loop
    for dummy=0:1:nber_radii-1    %rotate until we have rotated by 2pi
        if flag==1
            break
        else
            var_radius=[1 0];       %unit vector aligned with the x-axis (positve)
            Rot_matrix= [cos(angle_incr*dummy) -sin(angle_incr*dummy);...
                sin(angle_incr*dummy) cos(angle_incr*dummy)];
            
            var_radius = var_radius * Rot_matrix ;
            var_radius_inc = var_radius * Rot_matrix;
            %Check value of pixel at end of projected ray
            pixel_coord=round(c_2D+var_radius);
            
            %Store the size of var_radius when edge is reached (this can be improved)
            while slice2D_temp(pixel_coord(1), pixel_coord(2),row)>45
                if max(pixel_coord)==slice2D_size*2 || min(pixel_coord)==1
                    dummy2=dummy*angle_incr;
                    radius_error=sprintf...
                        ('edge is out of bound at %d th centrepoint and %d degrees ',row,dummy2);
                    disp(radius_error)
                    stored_radii(dummy+1)=norm(var_radius);
                    %dummy=nber_radii-1;
                    %flag=1;
                    pause
                    %var_radius=0; %set to zero to identify error
                    break
                else
                    var_radius= var_radius+var_radius_inc;
                    pixel_coord=round(c_2D+var_radius);
                end
            end
            stored_radii(dummy+1)=norm(var_radius);
        end
        
    end
    %Store the measure as stats in a row of the higher stroing matrix
    Stored_radii_para= cat(1,Stored_radii_para,...
        [mean(stored_radii), max(stored_radii),min(stored_radii)]);
    
end
figure (10)
row=1:1:size(F3,1);
plot(row, Stored_radii_para(:,1),'r')    %plot mean
hold on
plot(row, Stored_radii_para(:,2),'g')    %plot max
plot(row, Stored_radii_para(:,3))        %plot min in blue
%plot(dummyk,(10+cos(dummyk/6)),'k') %plot expected value in black


