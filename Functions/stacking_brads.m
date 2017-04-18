% stacking 2 of .lsm sequence
clear all
close all
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Capillaries work\Images')
%% First Stack - magenta
% Load and import all images into one 3D matrix
Stack = [];     % Set empty matrix
%t=50/256;    % Set threshold
for i=1:1:9
    Stack = cat(3,Stack, imread(sprintf('26_z00%d_c001.tif',i)));
end

for i=10:1:42
   Stack = cat(3,Stack, imread(sprintf('26_z0%d_c001.tif',i)));
end

Stack_p= Stack(601:700,601:700,5:38);

x=linspace (1,size(Stack_p,1),size(Stack_p,1));
y=linspace(1,size(Stack_p,2),size(Stack_p,2));
z=linspace (1,size(Stack_p,3),size(Stack_p,3));

figure
fv=isosurface(x,y,z,Stack_p,30);
p = patch(fv); 
%isonormals(M,p)
p.FaceColor = 'magenta';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight 
lighting gouraud
drawnow

%input dbcont to go to next step
%keyboard

%% Second Stack - red
Stack2 = [];     % Set empty matrix
%t=50/256;    % Set threshold
for i=1:1:9
    Stack2 = cat(3,Stack2, imread(sprintf('26_z00%d_c002.tif',i)));
end

for i=10:1:42
   Stack2 = cat(3,Stack2, imread(sprintf('26_z0%d_c002.tif',i)));
end

Stack2_p= Stack2(601:700,601:700,5:38);

hold on
fv=isosurface(x,y,z,Stack2_p,30);
p = patch(fv); 
%isonormals(M,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight 
lighting gouraud
drawnow

%input dbcont to go to next step
%keyboard

%% Third Stack - 
Stack3 = [];     % Set empty matrix
%t=50/256;    % Set threshold
for i=1:1:9
    Stack3 = cat(3,Stack3, imread(sprintf('26_z00%d_c003.tif',i)));
end

for i=10:1:42
   Stack3 = cat(3,Stack3, imread(sprintf('26_z0%d_c003.tif',i)));
end

Stack3_p= Stack3(601:700,601:700,5:38);

hold on
fv=isosurface(x,y,z,Stack3_p,30);
p = patch(fv); 
%isonormals(M,p)
p.FaceColor = 'green';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight 
lighting gouraud

%input dbcont to go to next step
%keyboard

%% Fourth Stack - 
%Stack4 = [];     % Set empty matrix
%t=50/256;    % Set threshold
%for i=1:1:9
%   Stack4 = cat(3,Stack4, imread(sprintf('26_z00%d_c004.tif',i)));
%end

%for i=10:1:42
%   Stack4 = cat(3,Stack4, imread(sprintf('26_z0%d_c004.tif',i)));
%end

%Stack4_p= Stack4(1:150,1:150,1:42);

%hold on
%fv=isosurface(x,y,z,Stack4_p,30);
%p = patch(fv); 
%isonormals(M,p)
%p.FaceColor = 'green';
%p.EdgeColor = 'none';
%daspect([1,1,1])
%/.view(3); axis tight
%camlight 
%lighting gouraud

