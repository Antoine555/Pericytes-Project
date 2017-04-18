xyradius=zeros((size(F3,1)-4),2); % matrix for recording the two point
    v= [1,0,0];%initialise direction for ray casting
    diameter=zeros((size(F3,1)-4),1);
    % Find the radii for each centrepopint
    fprintf([ 'Total number of centrepoints: ' num2str(size(F3,1)-4) ...
        '. Measurements taken:  '])
    for row=3:1:size(F3)-2
        
        if (row-3)<=10
            fprintf('\b%d',(row-3)); % update counter display
        else 
            fprintf('\b\b%d',(row-3)); % update counter display
        end

        % Find direction vector using derivative at point
        d=dF3dt(row,:);
        d=d/max(d)*50;
        if d(1)==0 && d(2)==0 % cross section is on a z-plane
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
        xyradius(row-2,1)=norm(vr3-F3(row,:)); % record the corresponding rad
        
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
        xyradius(row-2,2)=norm(vr3-F3(row,:)); % record the corresponding rad
        diameter(row-2)= xyradius(row-2,1)+xyradius(row-2,2);
    end
    fprintf('\n')
    disp( 'Finished measuring.')