%% Inputs: 
%             WhiteImg : path of an image taken through a diffuser, or of an entirely white scene
% 
%     GridModelOptions : struct controlling the model-bulding options
%           .ApproxLensletSpacing : A rough initial estimate to initialize the lenslet spacing
%                                   computation
%           .FilterDiskRadiusMult : Filter disk radius for prefiltering white image for locating
%                                   lenslets; expressed relative to lenslet spacing; e.g. a
%                                   value of 1/3 means a disk filte with a radius of 1/3 the
%                                   lenslet spacing
%                        .CropAmt : Image edge pixels to ignore when finding the grid
%                       .SkipStep : As a speed optimization, not all lenslet centers contribute
%                                   to the grid estimate; <SkipStep> pixels are skipped between
%                                   lenslet centers that get used; a value of 1 means use all
%           [optional] .Precision : 'single' or 'double'
% 
%  [optional] DebugDisplay : enables a debugging display, default false
function [LensletGridModel, GridCoords]=FPCUtilProcessWhiteImage(WhiteImageFullPath)

%% Load White image-----------------------------------------------------------------------------
WhiteImg=imread(WhiteImageFullPath);

%% initialize GridModelOption
GridModelOptions.ApproxLensletSpacing =23.316 ;
GridModelOptions = LFDefaultField( 'GridModelOptions', 'Orientation', 'horz' );
GridModelOptions = LFDefaultField( 'GridModelOptions', 'FilterDiskRadiusMult', 1/3 );
GridModelOptions = LFDefaultField( 'GridModelOptions', 'CropAmt', 25 );
GridModelOptions = LFDefaultField( 'GridModelOptions', 'SkipStep', 1 );

%% Display Hotpixel on the WhiteImage
DebugDisplay=false;

%% Build Hot Pixel Grid Coordinates
[LensletGridModel, GridCoords] = LFBuildLensletGridModel( WhiteImg, GridModelOptions, DebugDisplay );
end