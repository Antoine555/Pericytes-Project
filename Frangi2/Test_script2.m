%Frangi filter script

load('M2.mat');

%apply a threshold and close small regions
M3=M2>172;
M3=bwareaopen(M3,50,26);

% Frangi Filter the volume M3
options.BlackWhite=false;
options.FrangiScaleRange=[1 8];
M3filtered=FrangiFilter3D(M3,options);

% scale the intensities of the Frangi output
M3filt=M3filtered.*(10^9);
M3filt=uint8(M3filt); %transform to a uint8 type

%set figure and axes
figure
x=linspace (1,size(M2,1),size(M2,1));
y=linspace (1,size(M2,2),size(M2,2));
z=linspace (1,size(M2,3),size(M2,3));

%plot the volume M2 first, showing surface of intensity 0.5 in red
fv=isosurface(x,y,z,M3filt,20);
p = patch(fv);
%isonormals(M,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud

%overlay the new filtered volume with surface intensity i in blue
i=10;
fv=isosurface(x,y,z,M2filt,i);
p = patch(fv);
%isonormals(M2filtered,p)
p.FaceColor = 'blue';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud

   