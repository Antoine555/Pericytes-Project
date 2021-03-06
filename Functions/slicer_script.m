%slicer script
% Make sure the stack is a bouble array
%Stack3_pdb=double(Stack3_p(1:20,1:20,1:20));

%D=double(D);
xmin = min(x(:)); 
ymin = min(y(:)); 
zmin = min(z(:));

xmax = max(x(:)); 
ymax = max(y(:)); 
zmax = max(z(:));

%Create a plane
figure
hslice = surf(linspace(xmin,xmax,20),linspace(ymin,ymax,20),zeros(20));

rotate(hslice,[-1,0,0],-45);
xd = get(hslice,'XData');
yd = get(hslice,'YData');
zd = get(hslice,'ZData');

delete(hslice)

figure
colormap(jet)
h = slice(x,y,z,D,xd,yd,zd);
h.FaceColor = 'interp';
h.EdgeColor = 'none';
h.DiffuseStrength = 0.8;

%store the values in h
xh=h.XData;
yh=h.YData;
zh=h.ZData;
ch=h.CData;

hold on
hx = slice(x,y,z,Stack3_pdb,xmax,[],[]);
hx.FaceColor = 'interp';
hx.EdgeColor = 'none';

hy = slice(x,y,z,Stack3_pdb,[],ymax,[]);
hy.FaceColor = 'interp';
hy.EdgeColor = 'none';

hz = slice(x,y,z,Stack3_pdb,[],[],zmin);
hz.FaceColor = 'interp';
hz.EdgeColor = 'none';

daspect([1,1,1])
axis ([0, size(Stack3_pdb,1), 0, size(Stack3_pdb,2), 0, size(Stack3_pdb,3)])
view(-38.5,16)
camzoom(1.4)
camproj perspective