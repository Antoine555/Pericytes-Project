%% Recalculate centroid for 2D cross sections
    seg_slice2D_temp=zeros(41,41,size(F3,1));
    c_2D=zeros(size(F3,1),2);
    for row=1:size(F3,1)
        seg_slice2D_temp(:,:,row)= slice2D_temp(:,:,row)>150;
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
        I_temp=slice2D_temp(:,:,row);
        idx = s.PixelIdxList;
        sum_region = sum(I_temp(idx)-170*ones(size(idx,1),1));
        x_region = s.PixelList(:, 1);
        y_region = s.PixelList(:, 2);
        
        c_2D(row,1) = sum(x_region .* (double(I_temp(idx))...
            -170*ones(size(idx,1),1))) / sum_region;
        c_2D(row,2) = sum(y_region .* (double(I_temp(idx))...
            -170*ones(size(idx,1),1))) / sum_region;
        
    end
    
    figure(9)
    for row=1:size(F3,1)
        subplot(round((size(F3,1))/8+1),8,row)
        imagesc(seg_slice2D_temp(:,:,row))
    end
    
    pause