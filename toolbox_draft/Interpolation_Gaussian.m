% Image is interpolated to a defined scale (>1) from 381x383 original
% light-field dataset
load('Cocktails__Decoded.mat', 'LF');
LF=double(LF);        % convert LF to double
LF= LF/max(LF(:));

% Image is interpolated to a defined scale (>1) from 381x383 original light-field dataset
[NumberofRowSubImage,NumberofColumnSubImage, NumberOfMicroLensY, NumberOfMicroLensX, color]=size(LF);
scale=3;
interpolated_image=zeros(NumberOfMicroLensY*scale,NumberOfMicroLensX*scale);
sigma=.5;
%all neighbor pixels within a box R=half_windows_size are taken into considered
half_windows_size=2; 
% get central view  
Mid_NumberofRowSubImage=round(NumberofRowSubImage/2);
Mid_NumberofColumnSubImage=round(NumberofColumnSubImage/2);
tic;
for row=100:501         %2*scale:NumberOfRow*scale-2*scale
    for column=100:501  %2*scale:NumberOfRow*scale-2*scale
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
                    weigh_sum=weigh_sum+1/(2*pi()*sigma^2)*exp(-(distance.^2)/(2*sigma^2));
                    %if they are inside the boxs, sum up all cross-product GaussianWeigh*pixel_value         
                    temp=temp +1/(2*pi()*sigma^2)*exp(-(distance.^2)/(2*sigma^2))*...
                        LF(Mid_NumberofRowSubImage+pixel_r,Mid_NumberofColumnSubImage+pixel_c, r+pixel_r, c+pixel_c, 1);
                elseif or ( abs(pixel_r-interpolated_pixel_r)==half_windows_size,...
                        abs(pixel_c-interpolated_pixel_c)==half_windows_size)                
                    %calculate  distance
                    distance=sqrt((pixel_r-interpolated_pixel_r).^2+(pixel_c-interpolated_pixel_c).^2);
                    %sum weigh
                    weigh_sum=weigh_sum+1/(2*pi()*sigma^2)*exp(-(distance.^2)/(2*sigma^2));
                    %if they are inside the boxs, sum up all cross-product GaussianWeigh*pixel_value         
                    temp=temp +1/(2*pi()*sigma^2)*exp(-(distance.^2)/(2*sigma^2))*...
                        LF(Mid_NumberofRowSubImage+pixel_r,Mid_NumberofColumnSubImage+pixel_c, r+pixel_r, c+pixel_c, 1);
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