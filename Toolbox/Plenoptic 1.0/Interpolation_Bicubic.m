%Bicubic interpolation

load('C:\LFToolbox0.3_Samples\Images\Illum\IMG_0055__Decoded.mat', 'LF');
% calibrate LF 
LF=double(LF);        
LF= LF/max(LF(:));
tic;
% Image is interpolated to a defined scale (>1) from X x Y original light-field dataset
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
%initiate kernel variable
        weigh_sum=0;
        temp=0;
        %scale-down COORDINATE of a pixel in interpolated_image back to LF
        %size
        r=round(row/scale); 
        c=round(column/scale);
        %equivalently its COORDINATE inside u x v subimage is (interpolated_pixel_c,interpolated_pixel_r)
        interpolated_pixel_r= (row/scale-round(row/scale));
        interpolated_pixel_c= (column/scale-round(column/scale));
        %scan all pixels of subImage, which are inside windows
        for pixel_r=-half_windows_size:half_windows_size
            for pixel_c=-half_windows_size:half_windows_size
                if or ( abs(pixel_r-interpolated_pixel_r)>half_windows_size,...
                        abs(pixel_c-interpolated_pixel_c)>half_windows_size)                
                    %calculate  distance
                    distance=sqrt((pixel_r-interpolated_pixel_r).^2+(pixel_c-interpolated_pixel_c).^2);
                    %sum weigh
                    weigh_sum=weigh_sum+CubicConvolutionKernel(distance);
                    temp=temp+CubicConvolutionKernel(distance)*...
                        LF(Mid_NumberofRowSubImage+pixel_r,Mid_NumberofColumnSubImage+pixel_c, r+pixel_r, c+pixel_c, 1);;
                elseif or ( abs(pixel_r-interpolated_pixel_r)==half_windows_size,...
                        abs(pixel_c-interpolated_pixel_c)==half_windows_size)                
                    %calculate  distance
                    distance=sqrt((pixel_r-interpolated_pixel_r).^2+(pixel_c-interpolated_pixel_c).^2);
                    %sum weigh
                    weigh_sum=weigh_sum+CubicConvolutionKernel(distance);
                    temp=temp+CubicConvolutionKernel(distance)*...
                        LF(Mid_NumberofRowSubImage+pixel_r,Mid_NumberofColumnSubImage+pixel_c, r+pixel_r, c+pixel_c, 1);;
                else
                    continue
                end
            end
        end
        temp=temp/weigh_sum;
        interpolated_image(row,column)=temp;
    end
end
toc;
%print out interpolated image
interpolated_image=interpolated_image/max(interpolated_image(:));
imshow(interpolated_image,[])
imwrite(interpolated_image,'D:\Thesis_ light field_D\my works\Matlabcode\Results\interpolated_bicubic.jpg');