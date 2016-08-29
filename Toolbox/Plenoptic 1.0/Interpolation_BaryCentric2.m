%Barycentric interpolation
% interpolating_point(X) is calculated base on 3 neighbor_pixels, which
% form a matrix.
%
% (x1,y1)        X  
%
%                 o
% (x3,y3)  X            X  (x3,y3)
%
% Value of "x" point: g(x,y)=?1*f(x1,y1)+?2*f(x2,y2)+?3*f(x3,y3);
% where ?1,?2,?3 should satisfied system of equations: 
%        x = x1*?1+ x2*?2+ x3*?3
%        y = y1*?1+ y2*?2+ y3*?3
%        1 =    ?1+    ?2+    ?3

load('C:\LFToolbox0.3_Samples\Images\Illum\IMG_0055__Decoded.mat', 'LF');
% calibrate LF 
LF=double(LF);        
LF= LF/max(LF(:));
tic;
% Image is interpolated to a defined scale (>1) from 381x383 original light-field dataset
[NumberofRowSubImage,NumberofColumnSubImage, NumberOfMicroLensY, NumberOfMicroLensX, color]=size(LF);
scale=3;
interpolated_image=zeros(NumberOfMicroLensY*scale,NumberOfMicroLensX*scale);
%all neighbor pixels within a box R=half_windows_size are taken into considered
half_windows_size=1;    
% get central view  
Mid_NumberofRowSubImage=round(NumberofRowSubImage/2);
Mid_NumberofColumnSubImage=round(NumberofColumnSubImage/2);
tic;
for row=400:875         %2*scale:NumberOfMicroLensY*scale-2*scale            %Row number of interpolated image
    for column=550:1475  %2*scale:NumberOfMicroLensX*scale-2*scale     %Column number of interpolated image
        lambda=BarycentricCoefficience( floor(column/scale),floor(row/scale),...
                                        floor(column/scale),floor(row/scale)+1,...
                                        floor(column/scale)+1, floor(row/scale)+1,...
                                             column/scale,       row/scale);
        temp=    lambda(1)*LF(Mid_NumberofRowSubImage-1,Mid_NumberofColumnSubImage-1,floor(row/scale),floor(column/scale),1)...
                +lambda(2)*LF(Mid_NumberofRowSubImage-1,Mid_NumberofColumnSubImage+1,floor(row/scale),floor(column/scale),1)...
                +lambda(3)*LF(Mid_NumberofRowSubImage+1,Mid_NumberofColumnSubImage+1,floor(row/scale),floor(column/scale),1);
        interpolated_image(row,column)=temp;
    end
end
toc;
%print out interpolated image
interpolated_image=interpolated_image/max(interpolated_image(:));
imshow(interpolated_image,[])
imwrite(interpolated_image,'D:\Thesis_ light field_D\my works\Matlabcode\Results\interpolated_Barycentric.jpg');