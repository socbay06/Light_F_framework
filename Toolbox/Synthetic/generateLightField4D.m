%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Load a series of images, and save them along with additional data for
%   light field processing in a .mat file.
%
%   Gordon Wetzstein [gordonw@media.mit.edu]
%
%   University of British Columbia | PSM Lab | May 2011
%   MIT Media Lab | Camera Culture | Jan 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clf;

% path to output light field and folder containing the scene
%outpath = 'dice';
outpath = 'dragon';
%outpath = 'butterfly';
%outpath = 'messerschmitt';
%outpath = 'lucy';
%outpath = 'mini';
%outpath = 'xyzrgb_dragon';
%outpath = 'happy_buddha';
%outpath = 'teapot';
%outpath = 't-rex';
%outpath = 'bunnies';
%outpath = 'fishi';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bIs7x7LightField = false;
if  (strcmp(outpath, 'dice')==1) || ...
    (strcmp(outpath, 'dragon')==1) || ...
    (strcmp(outpath, 'butterfly')==1) || ...
    (strcmp(outpath, 'messerschmitt')==1) || ...
    (strcmp(outpath, 'lucy')==1) || ...
    (strcmp(outpath, 'mini')==1) || ...
    (strcmp(outpath, 'xyzrgb_dragon')==1) || ...
    (strcmp(outpath, 'happy_buddha')==1)
    bIs7x7LightField = true;
end


% size of the light field [y x] [mm]
if bIs7x7LightField
    lightFieldSize          = [75 100];
    lightFieldViewerDist    = inf;
else
    lightFieldSize          = [296.1 473.76];
    lightFieldViewerDist    = 2370;
end

% resolution in [#anglesY @anglesX #pixelsY #pixelsX #colorchannels]
% images will be resized or converted to grayscale depending on these
% parameters
if bIs7x7LightField
    lightFieldResolution    = [7 7 384 512 3];
else
    lightFieldResolution    = [5 5 525 840 3];
end

% path to images
datapath                = ['../' outpath '/'];
if ~bIs7x7LightField
    datapath = [datapath '/FOV20/'];
end

% field of view in degrees [y x]
if bIs7x7LightField
    fov                     = [10 10];
else
    fov                     = [20 20];
end

% max angle in both ways in v = tan(theta)
lightFieldAnglesX = -tan(pi*fov(2)/360):2*tan(pi*fov(2)/360)/(lightFieldResolution(1)-1):tan(pi*fov(2)/360);
lightFieldAnglesY = -tan(pi*fov(1)/360):2*tan(pi*fov(1)/360)/(lightFieldResolution(2)-1):tan(pi*fov(1)/360);
                    
% origin of light field in world space [y x z], x & y are defined in the 
% center of the images 
lightFieldOrigin        = [0 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% basename for files with count >10 and <100
filename = [datapath outpath '-'];
if lightFieldResolution(1)*lightFieldResolution(2)>10
   filename = [filename '0'];
end
filetype = '.png';

% initialize light field
lightField = zeros(lightFieldResolution);

count = 1;
for ky=1:lightFieldResolution(1)
    for kx=1:lightFieldResolution(2)
                        
        % filename of current image
        currentFilename = [filename num2str(count) filetype];        
        if count > 9
            currentFilename = [filename(1:end-1) num2str(count) filetype];            
        end
        
        % load image
        I = im2double(imread(currentFilename));
        
        % convert to grayscale if desired
        if lightFieldResolution(5) == 1
            I = rgb2gray(I);
        end
               
        % sort into light field datastructure
        lightField(ky,kx,:,:,:) = reshape(I, [1 1 lightFieldResolution(3) lightFieldResolution(4) lightFieldResolution(5)]);
        
        % draw current image
        imshow(I); drawnow;

        % increment counter
        count = count + 1;
    end
end

img = drawLightField4D(lightField);
img = imresize(img, [lightFieldResolution(3) lightFieldResolution(4)], 'bicubic');

imshow(img);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([datapath 'LightField4D.mat'], 'lightField', 'lightFieldAnglesY', 'lightFieldAnglesX', 'lightFieldSize', 'lightFieldResolution', 'lightFieldOrigin', 'lightFieldViewerDist');
