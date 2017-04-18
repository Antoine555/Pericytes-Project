% All in one test tube Slicer
% High end script that creates a test tube, skeletonises and
% then finds the radius.
clear all
close all

%% Create Test Tube
D=zeros(70,70,70);
k=1;
r=10+1*cos(k/6);
m=[35,35];
for k=1:70
    r=10+1*cos(k/6);
    for i=1:70
        for j=1:70
            if norm([i,j]-m)<=r
                D(i,j,k)=1;
            else
                D(i,j,k)=0;
            end
        end
    end
end

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
view(3); axis ([0 100 0 100 0 100]);
camlight
lighting gouraud
%[x,y,z]=ind2sub(D)
%plot3(D(1),D(2),D(3))


%% Skeletonisation


%addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Fast Marching')
%skel=skeleton(seg_stack3);

addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Skeleton3D')
skel = Skeleton3D(D);

% Plotting
figure
col=[.7 .7 .8];
hiso = patch(isosurface(D,0),'FaceColor',col,'EdgeColor','none');
hiso2 = patch(isocaps(D,0),'FaceColor',col,'EdgeColor','none');
axis equal;%axis off;
lighting phong;
isonormals(D,hiso);
alpha(0.5);
set(gca,'DataAspectRatio',[1 1 1])
camlight;
hold on;

%for i=1:length(skel)
% L=skel{i};
%plot3(L(:,2),L(:,1),L(:,3),'-','Color',rand(1,3));
% end

w=size(skel,1);
l=size(skel,2);
h=size(skel,3);
[x,y,z]=ind2sub([w,l,h],find(skel(:)));
plot3(y,x,z,'.','Markersize',8,'MarkerFaceColor','b','Color','b');
set(gcf,'Color','White');
view(140,80)
drawnow

%% Cleaning and labelling the skeleton
% Now use Skel2Graph on the skeleton
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Skel2graph3D')

skel2=skel;

w = size(skel2,1);
l = size(skel2,2);
h = size(skel2,3);

% convert skeleton to graph structure
[A,node,link] = Skel2Graph3D(skel2,8);

% convert graph structure back to (cleaned) skeleton
skel3 = Graph2Skel3D(node,link,w,l,h);

% display result
figure;
hold on;
for i=1:length(node)
    x1 = node(i).comx;
    y1 = node(i).comy;
    z1 = node(i).comz;
    for j=1:length(node(i).links)    % draw all connections of each node
        if(node(i).conn(j)<1)
            col='b'; % branches are blue
        else
            col='r'; % links are red
        end;
        
        % draw edges as lines using voxel positions
        for k=1:length(link(node(i).links(j)).point)-1
            [x3,y3,z3]=ind2sub([w,l,h],link(node(i).links(j)).point(k));
            [x2,y2,z2]=ind2sub([w,l,h],link(node(i).links(j)).point(k+1));
            line([y3 y2],[x3 x2],[z3 z2],'Color',col,'LineWidth',3);
        end;
    end;
    
    % draw all nodes as yellow circles
    plot3(y1,x1,z1,'o','Markersize',9,...
        'MarkerFaceColor','y',...
        'Color','k');
end;
%axis image;axis off;
axis ([0 70 0 70 0 70])
set(gcf,'Color','white');
drawnow;
view(-17,46);

%Plot the cleaned skeleton
% Plotting
figure
col=[.7 .7 .8];
hiso = patch(isosurface(D,0),'FaceColor',col,'EdgeColor','none');
hiso2 = patch(isocaps(D,0),'FaceColor',col,'EdgeColor','none');
axis equal;axis off;
lighting phong;
isonormals(D,hiso);
alpha(0.5);
set(gca,'DataAspectRatio',[1 1 1])
camlight;
hold on;

%for i=1:length(skel3)
% L=skel3{i};
%plot3(L(:,2),L(:,1),L(:,3),'-','Color',rand(1,3));
% end

w=size(skel3,1);
l=size(skel3,2);
h=size(skel3,3);
[x,y,z]=ind2sub([w,l,h],find(skel3(:)));
plot3(y,x,z,'.','Markersize',8,'MarkerFaceColor','b','Color','b');
set(gcf,'Color','White');
view(140,80)
drawnow

%% Frangi filter for extracting the eigenvalues

addpath ('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Frangi2')

% Frangi Filter the stent volume
options.BlackWhite=false;
options.FrangiScaleRange=[10 20];
[Dfiltered,Scale,Vx,Vy,Vz]=FrangiFilter3D(D,options);

%   % Show maximum intensity plots of input and result
figure,
subplot(2,2,1), imshow(squeeze(max(D,[],2)),[])
subplot(2,2,2), imshow(squeeze(max(Dfiltered,[],2)),[])
subplot(2,2,3), imshow(D(:,:,50),[])
subplot(2,2,4), imshow(Dfiltered(:,:,50),[])

%% Radius measuring

addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Extract Slice')
%Extract slice at nth centrepoint
c=[x,y,z];
Stored_radii_para=[];%set matrix for storing mean max min radii at each centrepoint
figure
for row=15:60 %for all centrepoints
    % Find direction vector (more than one method)
    %from centrepoints only
    %d=c(row,:)-c(row-1,:);
    %from frangi eigenvectors at centrepoints
    d=[Vx(c(row,:)),Vy(c(row,:)),Vz(c(row,:))];
    
    %Take slice
    [slice2D,sliceInd,subX,subY,subZ]=extractSlice(D,c(row,1),c(row,2),c(row,3),d(1),d(2),d(3),20);
    %Plot slices as we circulate through the centrepoints
    surf(subX,subY,subZ,slice2D,'FaceColor','texturemap','EdgeColor','none');
    axis([1 size(D,1) 1 size(D,2) 1 size(D,3)]);
    drawnow;
    
    % Could have a few lines here to relocated the centrepoint
    
    %Finding the radii on each slice at a time
    a=0;                        %set angle for rotation
    nber_radii=4; %choose the number of radii to be measured per centrepoint
    angle_incr=360/nber_radii;
    stored_radii=zeros(1,nber_radii); %set an empty matrix for storing measured radii
    c_2D=(size(slice2D)-[1 1])/2+[1 1];         %set coordinates of the centrepoint
    flag=0; %set a flag to be able to break out of a nested while loop
    for dummy=0:1:nber_radii    %rotate until we have rotated by 2pi
        if flag==1;
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
            while slice2D(pixel_coord(1), pixel_coord(2))>0
                if max(pixel_coord)==size(slice2D,1) || min(pixel_coord)==1
                    radius_error=sprintf('edge is out of bound for %d th centrepoint',row);
                    disp(radius_error)
                    flag=1;
                    var_radius=0; %set radius to zero to identify error
                    stored_radii(dummy+1)=norm(var_radius);
                    break
                else
                    var_radius= var_radius+var_radius_inc;
                    pixel_coord=round(c_2D+var_radius);
                end
            end
        
        end
        stored_radii(dummy+1)=norm(var_radius);
    end
    %Store the measure as stats in a row of the higher stroing matrix
    Stored_radii_para= cat(1,Stored_radii_para,...
        [mean(stored_radii), max(stored_radii),min(stored_radii)]);
    
end

dummyk=16:1:61;
plot(dummyk, Stored_radii_para(:,1),'r')    %plot mean
hold on
plot(dummyk, Stored_radii_para(:,2),'g')    %plot max
plot(dummyk, Stored_radii_para(:,3))        %plot min in blue
plot(dummyk,(10+cos(dummyk/6)),'k') %plot expected value in black


