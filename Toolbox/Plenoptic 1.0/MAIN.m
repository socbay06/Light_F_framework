%% extract view

%% DEPTH map

%% refocus
LF_Decoded_Matrix_path='D:\Thesis_ light field_D\LFP examples\lytro 10\Images\raws\Beers__Decoded.mat';
alpha=1.5;
patchSize=1;
%refocused_image1= FocalStackRendering_shiftAndAdd4(LF_Decoded_Matrix_path, alpha);
refocused_image2= FocalStackRendering_Patchsize(LF_Decoded_Matrix_path, patchSize);
refocused_image3= FocalStackRendering_cp2(LF_Decoded_Matrix_path, alpha); %recommended alpha: 2 | 0.2  |
%% super-resolution