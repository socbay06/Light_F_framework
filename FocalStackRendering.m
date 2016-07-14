alpha=0.8;   
sigma=.5;
% convert LF to double
% divide LF to max(LF(:))
[NumberofRowSubImage,NumberofColumnSubImage, NumberOfMicroLensY, NumberOfMicroLensX, color]=size(LF);
refocus_image=zeros(NumberOfMicroLensY,NumberOfMicroLensX);

% get central view  
Mid_NumberofRowSubImage=round(NumberofRowSubImage/2);
Mid_NumberofColumnSubImage=round(NumberofColumnSubImage/2);
tic;
for X= 1 :NumberOfMicroLensX
    for Y= 1:NumberOfMicroLensY  
        
        for pixel_r=-3:3%-(Mid_NumberofRowSubImage-1):(Mid_NumberofRowSubImage-1)
            for pixel_c=-3:3%-(Mid_NumberofColumnSubImage-1):(Mid_NumberofColumnSubImage-1)
                
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
imshow(refocus_image,[]);
toc;