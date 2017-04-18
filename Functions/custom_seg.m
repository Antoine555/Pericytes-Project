function [ seg_image ] = custom_seg( Image,method,t)
%CUSTOM_SEG will segment a 3D image according to the input method
% to be run from the Brad images folder
% requires 

%inputs:    Image: 3D image to be segmented
%           method: can be 'thresh', 'graph, or 'level'
%           t: chosen threshold, connectivity or lambda

%output:    seg_image: segmented image

switch method
    case 'thresh'
        seg_image=Image>=t;
        
    case 'graph'
        addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Graph cut simple')
        seg_image = [];     % Set empty matrix
        for i=1:1:size(Image,3)
            seg=graphcuts(Image(:,:,i),t,max(max(max(Image))));
            seg_image = cat(3,seg_image,seg);
        end
        
    case 'level'
        addpath('C:\Users\Anto\Documents\Oxford Uni\4th year\4YP\Matlab functions\Level set ChanVese\sfm_chanvese_demo')
        
        %compile necessary files
        disp('Compiling sfm_chanvese_mex...');
        mex sfm_chanvese_mex.cpp llist.cpp energy3c.cpp lsops3c.cpp
        disp('Compilation complete.');
        disp('');
        
        img=M(1:25,1:25,1:25); %.*(4095/255);
        
        % load image
        %img = imread('airplane.png');
        
        % prepare initialization
        mask = zeros(size(Image));
        mask(2:size(Image,1),2:size(Image,2),2:size(Image,3)) = 1;
        
        % set input parameters
        lambda = t; %for example t=0.001;
        iterations = 400;
        
        % perform segmentation
        [seg_image] = sfm_chanvese(img,mask,iterations,lambda);
        
        
    otherwise
        disp('Input method not supported. Please choose one of: thresh, graph, level')
end

end

