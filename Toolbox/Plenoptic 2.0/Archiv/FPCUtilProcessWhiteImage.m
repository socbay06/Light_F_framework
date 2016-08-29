load('D:\Thesis_ light field_D\LFP examples\Raytrix\Demo_01_RawWhite.mat');
WhiteImg=Demo_01_RawWhite;
XSpacing= 23.0344;
YSpacing= 21.66;
CropAmt=0;
EstAngle= 0;
Orientation= 'horz';
%---------------------------------------------------------------------------------------------
LensletGridModel = struct('HSpacing',XSpacing, 'VSpacing',YSpacing, 'HOffset',CropAmt, ...
    'VOffset',CropAmt, 'Rot',-EstAngle, 'Orientation', Orientation, ...
    'FirstPosShiftRow', 0,'FirstPosShiftCol', 0 );
LensletGridModel.UMax = round( (size(WhiteImg,2)-CropAmt*2)/XSpacing );     %174
LensletGridModel.VMax = round( (size(WhiteImg,1)-CropAmt*2)/YSpacing);      %123
%---------------------------------------------------------------------------------------------
GridCoords=zeros(LensletGridModel.VMax,LensletGridModel.UMax,2);
%---------------------------------------------------------------------------------------------
if LensletGridModel.FirstPosShiftRow< LensletGridModel.HSpacing
            GridCoords(1,1,1)=LensletGridModel.FirstPosShiftRow;
            GridCoords(1,1,2)=LensletGridModel.FirstPosShiftCol;
end
%---------------------------------------------------------------------------------------------
for vv=1:LensletGridModel.VMax      %123 
    %-----------------------------------------------------------------------------------------
    %coord Y trong cot dau tien cua dong tiep theo
    if and(GridCoords(vv,1,1)<(LensletGridModel.HSpacing/2-1), vv<LensletGridModel.VMax)
        GridCoords(vv+1,1,1)=GridCoords(vv,1,1)+LensletGridModel.HSpacing/2;
    elseif vv<LensletGridModel.VMax
        GridCoords(vv+1,1,1)=LensletGridModel.FirstPosShiftRow;
    end
    %-----------------------------------------------------------------------------------------
    for uu=1:LensletGridModel.UMax      %174 
        GridCoords(vv+1,1,2)=GridCoords(vv,1,2)+LensletGridModel.VSpacing;    %coord Y dau tien cong tiep theo
        if uu<LensletGridModel.UMax 
            GridCoords(vv,uu+1,1)=GridCoords(vv,uu,1)+LensletGridModel.HSpacing;     %coord cua cac X trong cung 1 dong
            GridCoords(vv,uu,2)=GridCoords(vv,1,2);         %coord cua cac Y trong cung 1 dong            
        else    %dung cho cot cuoi cung
            GridCoords(vv,uu,1)=LensletGridModel.UMax*LensletGridModel.HSpacing;     %coord X cua cot cuoi cung
            GridCoords(vv,uu,2)=GridCoords(vv,1,2);         %coord cua cac Y trong cung 1 dong   
        end
    end
end
        

%--Find offset to nearest peak for each--
GridCoordsX = GridCoords(:,:,1);
GridCoordsY = GridCoords(:,:,2);

%BuildGridCoords = cat(3,GridCoordsX,GridCoordsY);
BuildGridCoords = round(GridCoords);
%xlswrite('GridCoords.xlsx',BuildGridCoords);
