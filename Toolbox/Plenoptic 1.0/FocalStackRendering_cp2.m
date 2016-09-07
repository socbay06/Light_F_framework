function refocused_image= FocalStackRendering_cp2(LF_Decoded_Matrix_path, alpha)
%% Usage:
% refocus image using geometrical model
% Input:
%       LF_Decoded_Matrix_path: location of the LF decoded matrix
%       alpha: refocus factor 
%           >1  refocus further point
%           <1  refocus a nearer point
% Output:
%       Refocus Image
%% add toolbox LF
LFMatlabPathSetup
%% get all metadata from .lfp file
load(LF_Decoded_Matrix_path);
LF=double(LF);        % convert LF to double
LF= LF/max(LF(:));    % divide LF to max(LF(:))
f=1000*WhiteImageMetadata.devices.lens.focalLength;  %millimeter
d=WhiteImageMetadata.devices.mla.lensPitch/WhiteImageMetadata.devices.sensor.pixelPitch;
F= 43; %millimeter
m=F+WhiteImageMetadata.devices.lens.infinityLambda; %distance between focal_plane_when_object_is_in_infinity and MLA
%% u v y x parameters
UMax= size(LF,2);
VMax= size(LF,1);
XMax= LensletGridModel.UMax;
YMax= LensletGridModel.VMax;
%% get Grid model
cp=LFBuildHexGrid(LensletGridModel);
%% calculate intermediate values
beta=F/(m-F);
r_beta= d*beta/2;
v_0=ceil(VMax/2);
u_0=ceil(UMax/2);
%% filtering
refocused_image=squeeze(LF(v_0,u_0,:,:,1));
tic;
for xx=3:XMax-3
    for yy=3:260   %3:YMax-3
    temp=0;
    i=0;
        for uu=-2:2
        for vv=-2:2
        for i= -2:2
        for j= -2:2                
            if abs( (1/alpha-1)*F/f*sqrt((vv).^2+(uu).^2)+beta*sqrt((cp(yy,xx,1)-cp(yy+j,xx+i,1)).^2+(cp(yy,xx,2)-cp(yy+j,xx+i,2)).^2) ) >2*r_beta
                continue
            end
            %else sum up
            temp=temp+LF(v_0+vv,u_0+uu,yy+j,xx+i,1);
            i=i+1;
        end
        end
        end
        end
        refocused_image(yy,xx)=temp/i;
    end
end 
toc;
%% output

%% display
figure;imshow(refocused_image,[])
