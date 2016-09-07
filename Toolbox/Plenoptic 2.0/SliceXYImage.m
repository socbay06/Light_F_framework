function LF= SliceXYImage(RawImageFullPath, LensletGridModel, GridCoords)
RawImg=imread(RawImageFullPath);
%Demo_01_Raw;

% SliceXYImage( LensletGridModel, LensletImage, WhiteImage, DecodeOptions )
XSize = LensletGridModel.UMax;
YSize = LensletGridModel.VMax;
MaxSpacing = (max(LensletGridModel.HSpacing, LensletGridModel.VSpacing));  % Enforce square output in s,t
USize = round(MaxSpacing) ; % force odd for centered middle pixel -- H,VSpacing are even, so +1 is odd
VSize = round(MaxSpacing) ;
halfSTSize=floor(MaxSpacing/2);
LF = zeros(VSize, USize, YSize, XSize,3);
GridCoords=round(GridCoords);

for yy=1:YSize      %1:131
    for xx=1:XSize  %1:172
        if ((xx*LensletGridModel.HSpacing<MaxSpacing) ||(yy*LensletGridModel.HSpacing<MaxSpacing)||...
            (xx*LensletGridModel.HSpacing==MaxSpacing) ||(yy*LensletGridModel.HSpacing==MaxSpacing)||...
            (size(RawImg,1)-yy*LensletGridModel.VSpacing<2)||...
            (size(RawImg,2)-xx*LensletGridModel.HSpacing<2) )
            continue
        end
        %else cut sub-image
        LF(:,:,yy,xx,1)=RawImg(GridCoords(yy,xx,2)+(-halfSTSize:halfSTSize),GridCoords(yy,xx,1)+(-halfSTSize:halfSTSize));
    end
end
end