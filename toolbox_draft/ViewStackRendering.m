
[NumberofRowSubImage,NumberofColumnSubImage, NumberOfMicroLensX, NumberOfMicroLensY, color]=size(LF);
view_x=1;
view_y=1;
viewMid_X=round(NumberofRowSubImage/2);
viewMid_Y=round(NumberofColumnSubImage/2);
view_out=[];
view_out=squeeze(LF((viewMid_X+view_x),(viewMid_Y+view_y), :, :, :));