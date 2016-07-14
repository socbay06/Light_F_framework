% Image is interpolated to a defined scale (>1) from 381x383 original
% light-field dataset
function Interpolation_Gaussian(LF)
interpolated_image=[];
[NumberofRowSubImage,NumberofColumnSubImage, NumberOfMicroLensX, NumberOfMicroLensY, color]=size(LF);
scale=5;
sigma=.5;
%all neighbor pixels within a box R=half_windows_size are taken into considered
half_windows_size=1;    

for row=2*scale:NumberOfMicroLensX*scale-2*scale
    for column=2*scale:NumberOfMicroLensY*scale-2*scale
        interpolated_image(row,column)=0;
        %scale-down COORDINATE of a pixel in interpolated_image back to LF
        %size
        r=round(row/scale); 
        c=round(column/scale);
        
        %equivalently its COORDINATE inside 9x9 subimage is (interpolated_pixel_c,interpolated_pixel_r)
        interpolated_pixel_r= (row/scale-round(row/scale))+5;
        interpolated_pixel_c= (column/scale-round(column/scale))+5;
        
        %scan all pixels of subImage, which are inside windows
        for pixel_r=5-half_windows_size:5+half_windows_size
            for pixel_c=5-half_windows_size:5+half_windows_size
                if or ( abs(pixel_r-interpolated_pixel_r)>half_windows_size,...
                        abs(pixel_c-interpolated_pixel_c)>half_windows_size)
                    continue
                end
        %if they are inside the boxs, sum up all cross-product GaussianWeigh*pixel_value         
                interpolated_image(row,column)=...
                    interpolated_image(row,column)...
                    +1/(2*pi()*sigma^2)*exp(-((pixel_r-interpolated_pixel_r).^2+(pixel_c-interpolated_pixel_c).^2)/(2*sigma^2))...
                    *LF(pixel_c,pixel_r, r+(5-pixel_r), c+(5-pixel_c), 1);
            end
        end
    end
end
%print out interpolated image
%imshow(interpolated_image,[])
end