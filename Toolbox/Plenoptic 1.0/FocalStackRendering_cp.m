function refocused_image= FocalStackRendering_cp(LF_Decoded_Matrix_path, alpha)
%% add toolbox LF
LFMatlabPathSetup
%% get all metadata from .lfp file
load(LF_Decoded_Matrix_path);
LF=double(LF);        % convert LF to double
LF= LF/max(LF(:));    % divide LF to max(LF(:))
f=WhiteImageMetadata.devices.lens.focalLength;
d=WhiteImageMetadata.devices.mla.lensPitch/WhiteImageMetadata.devices.sensor.pixelPitch;
F= 84; %millimeter
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
r_alphaBeta= alpha*r_beta;
v_0=ceil(VMax/2);
u_0=ceil(UMax/2);
%% filtering
refocused_image=zeros(YMax, XMax);
tic;        
parfor xx=5:XMax-5
    for yy=5:YMax-5
    %temp=zeros(YMax,XMax);
        for uu=-4:4
        for vv=-4:4
        for i= -4:4
        for j= -4:4                
            if abs( (alpha-1)*F/f*sqrt((vv).^2+(uu).^2)+alpha*beta*sqrt((cp(yy,xx,1)-cp(yy+j,xx+i,1)).^2+(cp(yy,xx,2)-cp(yy+j,xx+i,2)).^2) )==2*r_alphaBeta||...
               abs( (alpha-1)*F/f*sqrt((vv).^2+(uu).^2)+alpha*beta*sqrt((cp(yy,xx,1)-cp(yy+j,xx+i,1)).^2+(cp(yy,xx,2)-cp(yy+j,xx+i,2)).^2) )> 2*r_alphaBeta
                continue
            end
            %else sum up
            refocused_image(yy,xx)=LF(v_0+vv,u_0+uu,yy,xx,1);
        end
        end
        end
        end
        %refocused_image=cat(3,refocused_image,temp);
    end
end 
toc;
%% output
refocused_image=sum(refocused_image,3);
%refocused_image=refocused_image/(UMax*VMax);

%% display
%refocused_image(Y_0,X_0)=0;
figure;imshow(refocused_image,[])
SaveFileName=strcat('FocalStackRendering_cp_alpha_',num2str(alpha),'.png');
imwrite(refocused_image,SaveFileName)