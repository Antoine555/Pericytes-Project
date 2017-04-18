% Create a sphere of fixed (pixel size) 
s_centre=[26,26,26];
s_volume=zeros(51,51,51);
radius=3/0.415;
for i=1:51
    for j=1:51
        for k=1:51
            d= norm([i,j,k]-s_centre);
            if d<=radius
                s_volume(i,j,k)=1;
            else
            end             
        end
    end
end
