%% get all RAW data, normalize
load('D:\Thesis_ light field_D\LFP examples\lytro 10\Images\raws\Beers__Decoded.mat');
LF=double(LF);        % convert LF to double
LF= LF/max(LF(:));    % divide LF to max(LF(:))
%% u v s t parameters
UMax= size(LF,2);
VMax= size(LF,1);
XMax= size(LF,4);
YMax= size(LF,3);
v_0=ceil(VMax/2);
u_0=ceil(UMax/2);
%% refocus parameter
alpha=4;
refocused_image=[];
%% do something
refocused_image=zeros(YMax+2*ceil(VMax/alpha),XMax+2*ceil(UMax/alpha),UMax*VMax);
tic;
for uu=-3:3
    for vv=-3:3
        temp=squeeze(LF(v_0+vv,u_0+uu,:,:,1));
        temp=padarray(temp,[ceil(VMax/alpha) ceil(UMax/alpha)]);
        u_shift=round(uu/alpha);
        v_shift=round(vv/alpha);
        weigh=( (u_shift-uu/alpha).^2+(v_shift-vv/alpha).^2 );
        temp=1*circshift(temp,[v_shift u_shift]);
        refocused_image(:,:,(u_0+uu)*UMax+v_0+vv)=temp;
        temp=[];
    end
end
toc;
%% output
refocused_image=sum(refocused_image,3);
%% display
figure;imshow(refocused_image,[])
SaveFileName=strcat('FocalStackRendering_shiftAndAdd_alpha_',num2str(alpha),'.png');
imwrite(refocused_image,SaveFileName)
