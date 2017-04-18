function ExploreRGB(input, mode,opts)
% Explore RGB imagery, plane by plane, and as grayscale
%
% ExploreRGB facilitates the viewing of all color planes, and grayscale and
% RGB versions, of input RGB image. The different colorplanes or
% conversions can be "selected" for export to the base workspace once
% clicked to expand.
%
% exploreRGB(input)
%     Simplest syntax: Simply view grayscale and color-plane options.
%
% exploreRGB(input, mode)
%     Mode can be specified as either 'normal' or 'advanced.'
%     (Default is 'normal', and is equivalent to using the simpler syntax
%     above.)
%
%     'advanced' also shows visualization of the RGB image as HSV, YCbCr,
%     and LAB. In R2014b+, it uses the new interface rgb2lab, and also
%     includes XYZ colorspace.
%
%     Alternatively, 'mode' may be specified as:
%     1 (= 'normal'), or 3 (= 'advanced'). (Note: Mode 2 is deprecated; new
%     behavior of expandAxes obviates it.)
%
%
%%% EXAMPLE:
% rgbImg = imread('peppers.png');
% exploreRGB(rgbImg)
%
% exploreRGB(rgbImg,'advanced')
% exploreRGB(rgbImg,1)
%
% Written by Brett Shoelson, PhD
% brett.shoelson@mathworks.com
% 4/24/08
%
% Modified:
%   2/18/2010 to provide more 'mode' control. Implemented 'simple' mode,
%   and 'normal' as a default.
%   2/02/2011 to provide an output handle to the figure (for convenience of
%   closing the figure)
%   5/07/2012 fixed labeling on ycbcr (r/b -> b/r)
%   2/24/13 Provide default input for nargin == 0; allow for loading of
%           additional images from workspace or file; suppress figure menu
%   3/11/13 'simple' mode deprecated.
%   6/25/14 Set the figure's handlevisibility to 'callback' to avoid
%           unintentional overwrites. Ensured compatibility with R2014b.
%   9/26/14 Disabled writing of images to axes' appdata. (Not sure why I
%           was doing that; in retrospect, it seems unwise.)
%
% Copyright 2008 - 2014 The MathWorks, Inc.
%
% See Also: customGray, rgb2gray, rgb2hsv, rgb2ycbcr, makecform,
%           applycform, rgb2lab, rgb2xyz

if nargin < 1 || isempty(input)
	tmp = load('mandrill');
	input = ind2rgb(tmp.X,tmp.map);
end
if ischar(input)
	input = imread(input);
end
if size(input,3)~=3
	error('ExploreRGB: Input must be an RGB image.')
end

if nargin < 2
	mode = 'normal';
elseif isnumeric(mode)
	if mode == 1, mode = 'normal';
	elseif mode == 3, mode = 'advanced';
	else
		error('ExploreRGB: Unsupported mode.')
	end
end
if ~ismember(lower(mode),{'normal','simple','advanced'})
	error('ExploreRGB: If specified, ''mode'' must be one of: ''normal'',''simple'',''advanced''.')
end

singleton = true; %Change to false to allow multiple instances
cleanup = [];
if singleton
	if ~isempty(findall(0,'tag','exploreRGBFigure'))
		cleanup = getappdata(findall(0,'tag','exploreRGBFigure'),...
			'AppCleanupCode');
	end
	delete(findall(0,'tag','exploreRGBFigure'));
end
% Clear the current figure
exploreRGBFighandle = figure('numbertitle','off',...
	'name','ExploreRGB: Click axes to EXPAND; Shift-Click to GENERATE CODE',...
	'color',[0.1 0.1 0.1],...
	'windowstyle','normal',...
	'units','normalized',...
	'position',[0.1 0.1 0.8 0.8],...[0.2 0.2 0.6 0.6],...
	'menubar','none',...
	'tag','exploreRGBFigure');
if ~isempty(cleanup)
	setappdata(exploreRGBFighandle,...
		'AppCleanupCode', cleanup);
end
if nargin > 2
	% opts specfied
	set(exploreRGBFighandle,...
		'units',opts.Units,...
		'windowstyle',opts.WindowStyle,...
		'position',opts.Position);
end
figh = exploreRGBFighandle;
setappdata(exploreRGBFighandle,'mode',mode);
ht = uitoolbar(exploreRGBFighandle);
tmp = im2double(imread('file_open.png'));
tmp(tmp==0) = NaN;
uitoggletool(ht,...
	'CData', tmp,...
	'oncallback', @GetNewFile,...
	'offcallback', '',...
	'Tooltipstring', 'Load new image.',...
	'Tag', 'loadImageTool');
tmp = getIcon(mode);
tmt = uitoggletool(ht,...
	'CData', tmp,...
	'oncallback', @toggleMode,...
	'offcallback', '',...
	'Tooltipstring', 'Toggle ''mode'': simple->advanced/advanced->simple',...
	'Tag', 'toggleMode');
tmp = im2double(imread('C.png'));
uitoggletool(ht,...
	'CData', tmp,...
	'oncallback', @callCustomGray,...
	'offcallback', '',...
	'Tooltipstring', 'Create CUSTOM gray image',...
	'Tag', 'customGray');
% Subplots
if strcmp(mode,'advanced')
	if verLessThan('matlab','8.4')
		nrows = 5;
	else
		nrows = 6;
	end
else
	nrows = 2;
end
ax(1) = subplot(nrows,3,1.5);
imh(1) = imshow(input);
%setappdata(gca,'currIm',input);
title('Original (RGB)',...
	'color','w',...
	'fontsize',14);

ax(2) = subplot(nrows,3,2.5);
img = rgb2gray(input);
imh(2) = imshow(img);
%setappdata(gca,'currIm',img);
title('Grayscaled',...
	'color',[0.4 0.4 0.4],...
	'fontsize',14);

ind = 3;
ax(ind) = subplot(nrows,3,ind+1);
img = input(:,:,1);
imh(3) = imshow(img);
%setappdata(gca,'currIm',img);
title('Red Plane',...
	'color',[0.7 0 0],...
	'fontsize',14);

ind = ind+1;
ax(ind) = subplot(nrows,3,ind+1);
img = input(:,:,2);
imh(ind) = imshow(img);
%setappdata(gca,'currIm',img);
title('Green Plane',...
	'color',[0 0.7 0],...
	'fontsize',14);

ind = ind+1;
ax(ind) = subplot(nrows,3,ind+1);
img = input(:,:,3);
imh(ind) = imshow(img);
%setappdata(gca,'currIm',img);
title('Blue Plane',...
	'color','b',...
	'fontsize',14)

if strcmp(mode,'advanced')
	hsv = rgb2hsv(input);
	titleColor = [1 0 1]*0.6;
	ind = ind+1;
	ax(ind) = subplot(nrows,3,ind+1);
	img = hsv(:,:,1);
	imh(ind) = imshow(img);
	%setappdata(gca,'currIm',img);
	title('Hue',...
		'color',titleColor,...
		'fontsize',14)
	
	ind = ind+1;
	ax(ind) = subplot(nrows,3,ind+1);
	img = hsv(:,:,2);
	imh(ind) = imshow(img);
	%setappdata(gca,'currIm',img);
	title('Saturation',...
		'color',titleColor,...
		'fontsize',14)
	
	ind = ind+1;
	ax(ind) = subplot(nrows,3,ind+1);
	img = hsv(:,:,3);
	imh(ind) = imshow(img);
	%setappdata(gca,'currIm',img);
	title('Value',...
		'color',titleColor,...
		'fontsize',14)
	
	ycbcr = rgb2ycbcr(input);
	titleColor = [0 1 1]*0.6;
	ind = ind+1;
	ax(ind) = subplot(nrows,3,ind+1);
	img = ycbcr(:,:,1);
	imh(ind) = imshow(img);
	%setappdata(gca,'currIm',img);
	title('Luminance',...
		'color',titleColor,...
		'fontsize',14)
	
	ind = ind+1;
	ax(ind) = subplot(nrows,3,ind+1);
	img = ycbcr(:,:,2);
	imh(ind) = imshow(img);
	%setappdata(gca,'currIm',img);
	title('Blue Chrominance',...
		'color',titleColor,...
		'fontsize',14)
	
	ind = ind+1;
	ax(ind) = subplot(nrows,3,ind+1);
	img = ycbcr(:,:,3);
	imh(ind) = imshow(img);
	%setappdata(gca,'currIm',img);
	title('Red Chrominance',...
		'color',titleColor,...
		'fontsize',14)
	
	%LAB
	if verLessThan('matlab','8.4')
		cform = makecform('srgb2lab');
		lab = applycform(input,cform);
	else
		lab = rgb2lab(input);
		disp('Using new RGB2LAB interface; please read the doc for that function, as there are other options that may prove helpful.')
	end
	titleColor = [1 1 0]*0.6;
	ind = ind+1;
	ax(ind) = subplot(nrows,3,ind+1);
	img = lab(:,:,1);
	imh(ind) = imshow(img);
	%setappdata(gca,'currIm',img);
	title('L*a*b* -> L*',...
		'color',titleColor,...
		'fontsize',14)
	
	ind = ind+1;
	ax(ind) = subplot(nrows,3,ind+1);
	img = lab(:,:,2);
	imh(ind) = imshow(img);
	%setappdata(gca,'currIm',img);
	title('L*a*b* -> a*',...
		'color',titleColor,...
		'fontsize',14)
	
	ind = ind+1;
	ax(ind) = subplot(nrows,3,ind+1);
	img = lab(:,:,3);
	imh(ind) = imshow(img);
	%setappdata(gca,'currIm',img);
	title('L*a*b* -> b*',...
		'color',titleColor,...
		'fontsize',14)
	
	if ~verLessThan('matlab','8.4')
		disp('Showing results for new RGB2XYZ interface; please read the doc for that function, as there are other options that may prove helpful.')
		xyz = rgb2xyz(input);
		titleColor = [1 1 0]*0.6;
		ind = ind+1;
		ax(ind) = subplot(nrows,3,ind+1);
		img = xyz(:,:,1);
		imh(ind) = imshow(img);
		%setappdata(gca,'currIm',img);
		title('XYZ -> X',...
			'color',titleColor,...
			'fontsize',14)
		
		ind = ind+1;
		ax(ind) = subplot(nrows,3,ind+1);
		img = xyz(:,:,2);
		imh(ind) = imshow(img);
		%setappdata(gca,'currIm',img);
		title('XYZ -> Y',...
			'color',titleColor,...
			'fontsize',14)
		
		ind = ind+1;
		ax(ind) = subplot(nrows,3,ind+1);
		img = xyz(:,:,3);
		imh(ind) = imshow(img);
		%setappdata(gca,'currIm',img);
		title('XYZ -> Z',...
			'color',titleColor,...
			'fontsize',14)
	end
end
expandAxes(ax)
processStrings = {...
	'Original';
	sprintf('grayscale = rgb2gray(rgbImg);');
	sprintf('grayscale = rgbImg(:,:,1);');
	sprintf('grayscale = rgbImg(:,:,2);');
	sprintf('grayscale = rgbImg(:,:,3);');
	sprintf('grayscale = rgb2hsv(rgbImg);\ngrayscale = grayscale(:,:,1);');
	sprintf('grayscale = rgb2hsv(rgbImg);\ngrayscale = grayscale(:,:,2);');
	sprintf('grayscale = rgb2hsv(rgbImg);\ngrayscale = grayscale(:,:,3);');
	sprintf('grayscale = rgb2ycbcr(rgbImg);\ngrayscale = grayscale(:,:,1);');
	sprintf('grayscale = rgb2ycbcr(rgbImg);\ngrayscale = grayscale(:,:,2);');
	sprintf('grayscale = rgb2ycbcr(rgbImg);\ngrayscale = grayscale(:,:,3);');
	sprintf('cform = makecform(''srgb2lab'');\ngrayscale = applycform(rgbImg,cform);\ngrayscale = grayscale(:,:,1);');
	sprintf('cform = makecform(''srgb2lab'');\ngrayscale = applycform(rgbImg,cform);\ngrayscale = grayscale(:,:,2);');
	sprintf('cform = makecform(''srgb2lab'');\ngrayscale = applycform(rgbImg,cform);\ngrayscale = grayscale(:,:,3);');
	};
if ~verLessThan('matlab','8.4')
	processStrings{12} = sprintf('grayscale = rgb2lab(rgbImg);\n%%(Alternate: grayscale = rgb2lab(rgbImg,''WhitePoint'',''icc'');)\ngrayscale = grayscale(:,:,1);');
	processStrings{13} = sprintf('grayscale = rgb2lab(rgbImg);\n%%(Alternate: grayscale = rgb2lab(rgbImg,''WhitePoint'',''icc'');)\ngrayscale = grayscale(:,:,2);');
	processStrings{14} = sprintf('grayscale = rgb2lab(rgbImg);\n%%(Alternate: grayscale = rgb2lab(rgbImg,''WhitePoint'',''icc'');)\ngrayscale = grayscale(:,:,3);');
	processStrings{15} = sprintf('grayscale = rgb2xyz(rgbImg);\ngrayscale = grayscale(:,:,1);');
	processStrings{16} = sprintf('grayscale = rgb2xyz(rgbImg);\ngrayscale = grayscale(:,:,2);');
	processStrings{17} = sprintf('grayscale = rgb2xyz(rgbImg);\ngrayscale = grayscale(:,:,3);');
end
for ii = 1:numel(imh)
	iptaddcallback(imh(ii),...
		'buttondownfcn',{@dispProcess,processStrings{ii}});
end
shg
set(exploreRGBFighandle,...
	'handlevisibility','callback')
if nargout == 0
	clear figh
end

	function callCustomGray(varargin)
		if ~exist('customGray','file')
			beep
			fprintf('\nPlease download ''customGray'' from the File Exchange to use this feature:\n\nhttp://www.mathworks.com/matlabcentral/fileexchange/43716\n');
			return
		end
		customGray(input,exploreRGBFighandle)
	end

	function dispProcess(varargin)
		selType = get(ancestor(gcbo,'figure'),'selectiontype');
		if strcmp(selType,'extend')
			fprintf('\n%%*****\n%s\n%%*****\n',varargin{3});
		end
	end

	function GetNewFile(varargin)
		set(gcbo,'state','off');
		[imgin,cmap,fname,~,userCanceled] = getNewImage(false);
		if userCanceled
			return
		end
		if ~isempty(cmap)
			imgin = ind2rgb(imgin,cmap);
		end
		opts = get(exploreRGBFighandle);
		ExploreRGB(imgin,mode,opts);
	end

	function toggleMode(varargin)
		set(gcbo,'state','off');
		mode = getappdata(exploreRGBFighandle,'mode');
		tmp = getIcon(mode);
		set(tmt,'cdata',tmp);
		if strcmp(mode,'normal')
			mode = 'advanced';
		else
			mode = 'normal';
		end
		setappdata(exploreRGBFighandle,'mode',mode);
		opts = get(exploreRGBFighandle);
		ExploreRGB(input,mode,opts);
	end

	function icon = getIcon(mode,varargin)
		if strcmp(mode,'advanced')
			icon = [
				0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
				0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
				0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
				0 0 1 1 1 1 1 1 1 1 0 0 2 2 2 2 2 2 2 0 0;
				0 0 1 3 3 1 1 1 1 1 0 0 2 2 2 2 2 2 2 0 0;
				0 0 1 3 3 1 1 1 1 1 0 0 2 2 2 2 2 2 2 0 0;
				0 0 1 1 1 4 4 1 1 1 0 0 2 2 2 2 2 2 2 0 0;
				0 0 1 1 1 4 4 1 1 1 0 0 2 2 2 2 2 2 2 0 0;
				0 0 1 1 1 1 1 5 5 1 0 0 2 2 2 2 2 2 2 0 0;
				0 0 1 1 1 1 1 5 5 1 0 0 2 2 2 2 2 2 2 0 0;
				0 0 1 1 1 1 1 1 1 1 0 0 2 2 2 2 2 2 2 0 0;
				0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
				0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
				0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
				0 0 3 3 3 3 0 0 4 4 4 4 0 0 5 5 5 5 0 0 0;
				0 0 3 3 3 3 0 0 4 4 4 4 0 0 5 5 5 5 0 0 0;
				0 0 3 3 3 3 0 0 4 4 4 4 0 0 5 5 5 5 0 0 0;
				0 0 3 3 3 3 0 0 4 4 4 4 0 0 5 5 5 5 0 0 0;
				0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
				0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
				0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
			icon = label2rgb(icon,[1 1 1;0.6 0.6 0.6;1 0 0;0 1 0; 0 0 1],[0.15 0.15 0.15]);
		else
			icon = [
				0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
				0 0 0 3 3 3 3 3 3 0 0 0 2 2 2 2 2 2 0 0 0;
				0 0 0 4 4 4 4 4 4 0 0 0 2 2 2 2 2 2 0 0 0;
				0 0 0 5 5 5 5 5 5 0 0 0 2 2 2 2 2 2 0 0 0;
				0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
				0 0 0 3 3 3 0 0 0 4 4 4 0 0 0 5 5 5 0 0 0;
				0 0 0 3 3 3 0 0 0 4 4 4 0 0 0 5 5 5 0 0 0;
				0 0 0 3 3 3 0 0 0 4 4 4 0 0 0 5 5 5 0 0 0;
				0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
				0 0 0 6 6 6 0 0 0 6 6 6 0 0 0 6 6 6 0 0 0;
				0 0 0 6 6 6 0 0 0 6 6 6 0 0 0 6 6 6 0 0 0;
				0 0 0 6 6 6 0 0 0 6 6 6 0 0 0 6 6 6 0 0 0;
				0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
				0 0 0 7 7 7 0 0 0 7 7 7 0 0 0 7 7 7 0 0 0;
				0 0 0 7 7 7 0 0 0 7 7 7 0 0 0 7 7 7 0 0 0;
				0 0 0 7 7 7 0 0 0 7 7 7 0 0 0 7 7 7 0 0 0;
				0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
				0 0 0 8 8 8 0 0 0 8 8 8 0 0 0 8 8 8 0 0 0;
				0 0 0 8 8 8 0 0 0 8 8 8 0 0 0 8 8 8 0 0 0;
				0 0 0 8 8 8 0 0 0 8 8 8 0 0 0 8 8 8 0 0 0;
				0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
			icon = label2rgb(icon,[1 1 1;0.6 0.6 0.6;1 0 0;0 1 0; 0 0 1;1 0 1;0 1 1;1 1 0],[0.15 0.15 0.15]);
		end
	end

end
