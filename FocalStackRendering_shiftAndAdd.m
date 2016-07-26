alpha=2.05;   %1    %1.95   %
alpha=2.05;   %1    %1.95   %
sigma=.5;
load('C:\LFToolbox0.3_Samples\Images\F01\IMG_0009.RAW__Decoded.mat', 'LF');
LF=double(LF);        % convert LF to double
LF= LF/max(LF(:));    % divide LF to max(LF(:))
[NumberofRowSubImage,NumberofColumnSubImage, NumberOfMicroLensY, NumberOfMicroLensX, color]=size(LF);
refocus_image=zeros(NumberOfMicroLensY,NumberOfMicroLensX);

% get central view  
Mid_NumberofRowSubImage=round(NumberofRowSubImage/2);
Mid_NumberofColumnSubImage=round(NumberofColumnSubImage/2);
tic;
for X= 1 :NumberOfMicroLensX
    for Y= 1:NumberOfMicroLensY  
        
        for pixel_r=-2:2%-(Mid_NumberofRowSubImage-1):(Mid_NumberofRowSubImage-1)
            for pixel_c=-2:2%-(Mid_NumberofColumnSubImage-1):(Mid_NumberofColumnSubImage-1)
                outputImagePixelX=X+(Mid_NumberofRowSubImage+pixel_c)*(1-1/alpha);
                outputImagePixelY=Y+(Mid_NumberofColumnSubImage+pixel_r)*(1-1/alpha);
                if ((outputImagePixelY > NumberOfMicroLensY)||(outputImagePixelX > NumberOfMicroLensX)||...
                        (outputImagePixelX<1)||(outputImagePixelY<1))
                    continue
                end
                refocus_image(Y,X)=refocus_image(Y,X)+LF(Mid_NumberofColumnSubImage+pixel_c,...
                                             Mid_NumberofRowSubImage+pixel_r,... 
                                             round(outputImagePixelY),...
                                             round(outputImagePixelX),...
                                             1);
                end
            end
        end
end
refocus_image=refocus_image/max(refocus_image(:));
figure;
imshow(refocus_image,[]);
toc;