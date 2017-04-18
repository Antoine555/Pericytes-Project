%sequencing
figure
for i=1:36
    subplot(6,6,i), imshow(M2filtered(:,:,i),[]);
    i=i+1;
end

    