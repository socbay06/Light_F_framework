%take Cocktails__Decoded as example
load('Cocktails__Decoded.mat', 'LF');
load('Cocktails__Decoded.mat', 'LFMetadata');
LF=double(LF);        % convert LF to double
LF= LF/max(LF(:));
%input parameter:
ininfinityLambda=   LFMetadata.devices.lens.infinityLambda;
alpha=              2; %new plane to be refocus
focalLength=        1000*LFMetadata.devices.lens.focalLength;
ML_Diameter=        LFMetadata.devices.mla.lensPitch;
[NumberofRowSubImage,NumberofColumnSubImage, NumberOfMicroLensY, NumberOfMicroLensX, color]=size(LF);

%output image
refocus_image=zeros(NumberOfMicroLensY,NumberOfMicroLensX);

%start time counting
tic;

for X= 1 :NumberOfMicroLensX
    for Y= 1:NumberOfMicroLensY  
        
        for pixel_r=1:NumberofRowSubImage
            for pixel_c=1:NumberofColumnSubImage
                
            end
        end
    end
end

refocus_image=refocus_image/max(refocus_image(:));
figure;
imshow(refocus_image,[]);
toc;                
%---------------------------------
