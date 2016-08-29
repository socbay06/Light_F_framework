load('D:\Thesis_ light field_D\LFP examples\Raytrix\Demo_01_Raw.mat');
RawImg=Demo_01_Raw;

% SliceXYImage( LensletGridModel, LensletImage, WhiteImage, DecodeOptions )
USize = LensletGridModel.UMax;
VSize = LensletGridModel.VMax;
MaxSpacing = (max(LensletGridModel.HSpacing, LensletGridModel.VSpacing));  % Enforce square output in s,t
SSize = round(MaxSpacing) ; % force odd for centered middle pixel -- H,VSpacing are even, so +1 is odd
TSize = round(MaxSpacing) ;
halfSTSize=floor(MaxSpacing/2);
LF = zeros(TSize, SSize, VSize, USize,3);
GridCoords=round(GridCoords);

for vv=1:VSize      %1:123
    for uu=1:USize  %1:174
        if ((uu*LensletGridModel.HSpacing<MaxSpacing) ||(vv*LensletGridModel.HSpacing<MaxSpacing)||...
              (uu*LensletGridModel.HSpacing==MaxSpacing) ||(vv*LensletGridModel.HSpacing==MaxSpacing)||...
              (size(RawImg,1)-vv*LensletGridModel.HSpacing<2)||...
              (size(RawImg,2)-uu*LensletGridModel.HSpacing<MaxSpacing))
            continue
        end
        %else cut sub-image
        LF(:,:,vv,uu,1)=RawImg(GridCoords(vv,uu,2)+(-halfSTSize:halfSTSize),GridCoords(vv,uu,1)+(-halfSTSize:halfSTSize));
    end
end

IMGCheck=squeeze(LF(12,12,:,:,1));
%imagesc(IMGCheck)
%IMGCheck=RawImg(GridCoords(vv,uu,2)+(-halfSTSize:halfSTSize),GridCoords(vv,uu,1)+(-halfSTSize:halfSTSize));
figure; imshow(IMGCheck,[])
