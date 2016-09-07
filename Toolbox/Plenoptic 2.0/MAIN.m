LFMatlabPathSetup
%% decode whiteImage
WhiteImageFullPath='D:\Thesis_ light field_D\LFP examples\Raytrix\Demo_01_RawWhite.png';
[LensletGridModel, GridCoords]=FPCUtilProcessWhiteImage(WhiteImageFullPath);

%% store RAW image in 5D light field
RawImageFullPath='D:\Thesis_ light field_D\LFP examples\Raytrix\Demo_01_Raw.png';
LF=SliceXYImage(RawImageFullPath, LensletGridModel, GridCoords);

%% display yfor debugging
IMGCheck=squeeze(LF(12,12,:,:,1));
imshow(IMGCheck,[])
IMGFull = reshape(permute(LF,[1 3 2 4 5]),[size(LF,1)*size(LF,3) size(LF,2)*size(LF,4) 3]);
%imshow(IMGFull,[])