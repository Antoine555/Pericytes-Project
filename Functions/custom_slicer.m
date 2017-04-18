function [ s ] = custom_slicer(Stack,c,n, plane_size)
%CUSTOM_SLICER
%slices the image at the given centrepoint c and returns a 2D image
%direction/normal vector is found from previous centrepoint ((n-1)th)

%inputs
%Stack: original stack of double values
%c: centrepoint coordinates
%n: nth centrepoint
%plane_size: size (sides) of the cutting plane (set to bigger than the max
%diameter)

%x,y,z: cordinates of centrepoints c
%d: direction vector

%outputs:
%s: slice at c(n) of side plane_size and normal d

%Make sure the stack is an array of doubles
Stack=double(Stack);
x=c(:,1); y=c(:,2); z=c(:,3);

%Find size of stack
S=size(Stack);

xmin = 1;
ymin = 1;
zmin = 1;

xmax = S(1);
ymax = S(2);
zmax = S(3);

% Find direction vector (preceding i-th point)
%d=[x(i),y(i),z(i)]-[x(i-1),y(i-1),z(i-1)];     %find direction vector
d=c(n,:)-c(n-1,:);                              %find direction vector
%d=d/norm(d);                                   %normalise direction vector


%Create a points surface on the x-y plane
figure
hslice = surf(linspace(-plane_size/2,plane_size/2,plane_size),...
    linspace(-plane_size/2,plane_size/2,plane_size),zeros(plane_size));
hslice2=hslice;

%set exception points
pts1=0;
%pts2=0;

if d(1)==0 && d(2)==0 && d(3)==0
    pts1=pts1+1   %this exception should not happen/could be avoided in the upper script
%Avoid this in the upper script?    
%elseif dot(c)-c(i-1,:),c(i+1,:)-c(i,:))<= (1/sqrt(2))
%    pts2=pts2+1;    %points are not on the same edge of centreline
elseif d(3)==0 && d(2)==0 %plane on the z-y plane
    rotate(hslice2,[0,1,0],90);
elseif d(3)==0 && d(1)==0 %plane on the z-x plane
    rotate(hslice2,[1,0,0],90);
elseif d(3)==0; % plane does not cross z-axis
    beta=atan(d(1)/d(2)) *360/(2*pi);
    %Rx=[1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];
    rotate(hslice2,[1,0,0],90);
    rotate(hslice2,[0,0,1],beta);
    
elseif d(1)==0 && d(2)==0 %on the x-y plane
else %any other orientation of plane
%find angles from d, (x-axis (alpha) then z-axis (beta) rotations)
beta=atan(d(1)/d(2)) *360/(2*pi);
alpha=atan(d(3)/(d(1)^2+d(2)^2)) *360/(2*pi);

%create rotation matrices
%Rx=[1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];
%Rz=[cos(beta),-sin(beta),0;sin(beta),cos(beta),0;0,0,1];

%rotate so that d is the normal vector
rotate(hslice2,[1 0 0],-alpha);
rotate(hslice2, [0,0,1], -beta);

end

%translate plane to centrepoint
hslice2.XData=hslice2.XData+ones(size(hslice2.XData))*c(n,1);
hslice2.YData=hslice2.YData+ones(size(hslice2.YData))*c(n,2);
hslice2.ZData=hslice2.ZData+ones(size(hslice2.ZData))*c(n,3);

xd = get(hslice2,'XData');
yd = get(hslice2,'YData');
zd = get(hslice2,'ZData');

%delete(hslice)


%Plots
figure
colormap(jet)
h = slice(Stack,xd,yd,zd);
h.FaceColor = 'interp';
h.EdgeColor = 'none';
h.DiffuseStrength = 0.8;

s=h;

%store the values in h
xh=h.XData;
yh=h.YData;
zh=h.ZData;
ch=h.CData;

hold on
hx = slice(Stack,xmax,[],[]);
hx.FaceColor = 'interp';
hx.EdgeColor = 'none';


hy = slice(Stack,[],ymax,[]);
hy.FaceColor = 'interp';
hy.EdgeColor = 'none';

hz = slice(Stack,[],[],zmin);
hz.FaceColor = 'interp';
hz.EdgeColor = 'none';

daspect([1,1,1])
axis ([0, size(Stack,1), 0, size(Stack,2), 0, size(Stack,3)])
view(-38.5,16)
camzoom(1.4)
camproj perspective
end