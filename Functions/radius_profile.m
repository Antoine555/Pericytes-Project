
% measure radius by studying profile
Mprof=[];

% Measure the radii from these images
for row=3:1:size(F3,1)-4 %for all centrepoints
    nber_radii=4; %choose the number of radii to be measured per centrepoint
    angle_incr=360/nber_radii;
    stored_radii=zeros(1,nber_radii); %set an empty matrix for storing measured radii
    c_2D=slice2D_size*2/2+[1 1];   %set coordinates of the centrepoint
    
    for dummy=0:1:nber_radii-1    %rotate until we have rotated by 2pi
        
        
        var_radius=slice2D_size*[1 0];       %unit vector aligned with the x-axis (positve)
        Rot_matrix= [cos(angle_incr*dummy) -sin(angle_incr*dummy);...
            sin(angle_incr*dummy) cos(angle_incr*dummy)];
        
        var_radius = var_radius * Rot_matrix ;
        
        %Check value of pixel at end of projected ray
        pixel_coord=round(c_2D+var_radius);
        prof=improfile (slice2D_temp(:,:,row),[c_2D(1),pixel_coord(1)],...
            [c_2D(2),pixel_coord(2)] ,15);
        Mprof=[Mprof, prof];
        %Store the size of var_radius when edge is reached (this can be improved)
        
plot(Mprof(:,(row-3)*nber_radii+dummy+1),'r')    %plot mean
hold on

    end
    
       %plot min in blue

    
end
%Store the measure as stats in a row of the higher stroing matrix


