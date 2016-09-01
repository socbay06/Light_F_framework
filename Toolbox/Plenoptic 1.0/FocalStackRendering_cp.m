%% add toolbox LF
%LFMatlabPathSetup
%% get all metadata from .lfp file
load('IMG_0010.RAW__Decoded.mat');
LF=double(LF);        % convert LF to double
LF= LF/max(LF(:));    % divide LF to max(LF(:))
alpha=2;
f=WhiteImageMetadata.devices.lens.focalLength;
d=WhiteImageMetadata.devices.mla.lensPitch/WhiteImageMetadata.devices.sensor.pixelPitch;
F= 84; %millimeter
m=F+WhiteImageMetadata.devices.lens.infinityLambda; %distance between focal_plane_when_object_is_in_infinity and MLA
%% u v s t parameters
UMax= size(LF,2);
VMax= size(LF,1);
XMax= LensletGridModel.UMax;
YMax= LensletGridModel.VMax;
%% get Grid model
cp=LFBuildHexGrid(LensletGridModel);
%% calculate intermediate values
beta=m*F/(m-F);
r_beta= d*beta/2;
r_alphaBeta= alpha*r_beta;
v_0=ceil(VMax/2);
u_0=ceil(UMax/2);
%% choose focus pixel
X_0= 160;
Y_0= 217;
%% filtering
refocused_image=[];
tic;
for uu=u_0-4:u_0+4
    for vv=v_0-4:v_0+4
        temp=zeros(YMax,XMax);
        for xx=1:XMax
            for yy=1:YMax
                if abs( (1-1/alpha)*F/(f*beta*d)*sqrt((v_0-vv).^2+(u_0-uu).^2)+sqrt((cp(yy,xx,1)-cp(Y_0,X_0,1)).^2+(cp(yy,xx,2)-cp(Y_0,X_0,2)).^2) )==d||...
                   abs( (1-1/alpha)*F/(f*beta*d)*sqrt((v_0-vv).^2+(u_0-uu).^2)+sqrt((cp(yy,xx,1)-cp(Y_0,X_0,1)).^2+(cp(yy,xx,2)-cp(Y_0,X_0,2)).^2) )>d
                    continue
                end
                %else sum up
                temp(yy,xx)=LF(vv,uu,yy,xx,1);
            end
        end
        refocused_image=cat(3,refocused_image,temp);
    end
end 
toc;
%% output
refocused_image=sum(refocused_image,3);
%refocused_image=refocused_image/(UMax*VMax);

%% display
figure;imshow(refocused_image,[])
SaveFileName=strcat('FocalStackRendering_cp_alpha_',num2str(alpha),'.png');
%imwrite(refocused_image,SaveFileName)