function [img,cmap,fname,fpath,userCanceled] = getNewRGBImage(convertToDouble)
% Prompt for new RGB image from file or workspace
%
% SYNTAX:
% [img,map,fname,fpath,userCanceled] = getNewRGBImage(convertToDouble)
%
% INPUTS:
% convertToDouble: Binary (default = false)
% OUTPUTS:
% img, cmap: Image and corresponding colormap
% fname, fpath: filename, file path
% userCanceled: Logical, indicating whether the user canceled the operation
%
% Incorporates Scott Hirsch's UIGETVARIABLES; modified from getNewImage.
%
% Copyright MathWorks,Inc 2012
%
% Brett Shoelson, PhD
% brett.shoelson@mathworks.com
% 08/09/2013

% DEFAULTS
img  = [];cmap = [];fname = [];fpath = []; userCanceled = false;
if nargin == 0
    convertToDouble = false;
end
importopt = questdlg('Import from File or Workspace?','File Source','FILE','Workspace','Cancel','FILE');
switch importopt
    case 'Cancel'
        userCanceled = true;
        return
    case 'FILE'
        filterspec = imgformats(1);
        [fname,fpath] = uigetfile(filterspec, 'Select image file');
        if ~ischar(fpath)
            userCanceled = true;
            return
        end
        [img,cmap] = imread(fullfile(fpath,fname));
    case 'Workspace'
        prompts = 'MxNx3 Image';
        inputstring = 'Please select an RGB matrix representing your image.';
        % The convenience syntax for uigetvariables doesn't give enough
        % control, so specify a validator manually
        % * numeric
        % * 3D: third dimension must be length 3
        validateimage = @(im) isnumeric(im) &&((ndims(im)==3&&size(im,3)==3));%#ok
        img = uigetvariables(prompts,inputstring,validateimage);
        if isa(img,'cell')
            try
                img = img{1};
            end
        end
        if isempty(img)
            userCanceled = true;
        end
end
if ~userCanceled && convertToDouble
    img = im2double(img);
end
