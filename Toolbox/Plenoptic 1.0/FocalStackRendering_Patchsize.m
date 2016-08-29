alpha=2.05;   %1    %1.95   %
alpha=2.05;   %1    %1.95   %
sigma=.5;
load('C:\LFToolbox0.3_Samples\Images\F01\IMG_0009.RAW__Decoded.mat', 'LF');
LF=double(LF);        % convert LF to double
LF= LF/max(LF(:));    % divide LF to max(LF(:))
[NumberofRowSubImage,NumberofColumnSubImage, NumberOfMicroLensY, NumberOfMicroLensX, color]=size(LF);
refocus_image=[];
pixel_r=1;%-(Mid_NumberofRowSubImage-1):(Mid_NumberofRowSubImage-1)
pixel_c=1;%-(Mid_NumberofColumnSubImage-1):(Mid_NumberofColumnSubImage-1)

% get central view  
Mid_NumberofRowSubImage=round(NumberofRowSubImage/2);
Mid_NumberofColumnSubImage=round(NumberofColumnSubImage/2);
tic;
for Y= 1 :NumberOfMicroLensY
    temp2=[];
    for X= 1 :NumberOfMicroLensX  
        temp1=LF(   Mid_NumberofRowSubImage-pixel_r:Mid_NumberofRowSubImage+pixel_r,...
                    Mid_NumberofColumnSubImage-pixel_c:Mid_NumberofColumnSubImage+pixel_c,...
                    Y, X, 1);
        temp1=flip(temp1,2);
        temp1=flip(temp1,1);
        temp2=cat(2,temp2,temp1);                         
    end
    refocus_image=cat(1,refocus_image,temp2);
end
refocus_image=refocus_image/max(refocus_image(:));
figure;
imshow(refocus_image,[]);
toc;