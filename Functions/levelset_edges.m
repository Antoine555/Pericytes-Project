addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Level set toolbox\AOSLevelsetSegmentationToolboxM\data')
addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Level set toolbox\AOSLevelsetSegmentationToolboxM')
contour_ptM=[];
contour_pt_s=struct('cross_section',0,'coord',[]);
for row=3:1:size(F3-2)
I =  slice2D_temp(:,:,row);
phi = ac_SDF_2D('circle', size(I), c_2D(15,:),slice2D_size) ;
   %{   
smooth_weight = 1; image_weight = 10*1e-3; 
      delta_t = 2; n_iters = 1000; show_result = 1; 
      phi = ac_ChanVese_model(double(I), phi, smooth_weight, image_weight, ...
          delta_t, n_iters, show_result);
%}

     g = ac_gradient_map(I, 5); 
     mu = 100; propagation_weight = .05; GAC_weight = 5000; 
     delta_t = 1; n_ite = 100; show_result = 1; n_iters=100;
     phi = ac_hybrid_model(I-mu, phi, propagation_weight, GAC_weight, g, ...
       delta_t, n_iters, show_result);
close all
   c = contours(phi,[0,0]);
   contour_pt=[];
   idx = 1;
   while idx < size(c,2)
       n = c(2,idx);
   contour_pt = [contour_pt, c(:,idx+1:idx+n)];
   idx = idx+n+1;
   end
   contour_pt_s(row-2).cross_section=row;
   contour_pt_s(row-2).coord=permute(contour_pt,[2 1]);
end

   
   %% Connect edge points to form a spline perimeter
perimeter=struct('cross_section',0,'coord',zeros(size(contour_pt)));
figure (13)
 
for row=3:size(F3,1)-2
    
    perimeter(row-2).cross_section=row;
    perimeter(row-2).coord=[contour_pt_s(row-2).coord(:,:)];
    
    subplot(round((size(F3,1))/8+1),8,row)
    imagesc(slice2D_temp(:,:,row))
    hold on
    plot(contour_pt_s(row-2).coord(:,1),contour_pt_s(row-2).coord(:,2),...
        'Marker','.','LineStyle','none','MarkerEdgeColor','r','MarkerFaceColor','r')
    plot(c_2D(row,1),c_2D(row,2),'*')
end

pause

%% Calculate the enclosed surface

Areas=zeros(size(F3,1)-4,1);
for row=3:size(F3,1)-2
    Areas(row-2)=polyarea(perimeter(row-2).coord(:,1),...
        perimeter(row-2).coord(:,2));
    
end

%% Calculate the perimeter
Ps=zeros(size(F3,1)-4,1);
perimeteri=struct('cross_section',0,'coord',zeros(100,2));
for row=3:size(F3,1)-2
    perimeteri(row-2).cross_section=row;
    perimeteri(row-2).coord=interparc(100, perimeter(row-2).coord(:,1),...
        perimeter(row-2).coord(:,2));
    Ps (row-2)= 100*norm( ...
        [perimeteri(row-2).coord(3,1)-perimeteri(row-2).coord(2,1),...
        perimeteri(row-2).coord(3,2)-perimeteri(row-2).coord(2,2)]);
end

%% Plot interpolated perimeter
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

%% Recalculate the radius
Stored_radii_para2=[];
stored_vectors2=zeros(100,2);
for row=3:1:size(F3)-2
    stored_vetors2=perimeteri(row-2).coord-...
        [c_2D(row-2,1)*ones(100,1),c_2D(row-2,2)*ones(100,1)];
    stored_radii2=zeros(size(stored_vectors2,1),1);
    for point=1:1:100
        stored_radii2(point)=norm(stored_vetors2(point,:));
    end
    
    Stored_radii_para2= cat(1,Stored_radii_para2,...
        [mean(stored_radii2), max(stored_radii2),min(stored_radii2)]);
end

%% Plot radius, surface and perimeter along arclength
figure (15)
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