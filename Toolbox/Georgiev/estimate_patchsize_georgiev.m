% Working with focused dataset using Georgiev method
% Georgiev, T., & Lumsdaine, A. (2010). Focused plenoptic camera and rendering. Journal of Electronic Imaging, 19(2), 021106. http://doi.org/10.1117/1.3442712

lf=hdf5read('focused.flower.h5','/lf/focused/flower');
lf=double(lf);


[T,S,H,W,C] = size(lf);
% define search window
m = 31;

hm = floor((m-1)/2);
centerx = floor((W-1)/2)+1;
centery = floor((H-1)/2)+1;

% range of kx ky
kx = 1+hm:1:W-hm;
ky = 1+hm:1:H-hm;

% calculate the cost for left, right, up and down mimg seperately.
costxleft = zeros(T,S,size(kx,2))-1;
costxright = zeros(T,S,size(kx,2))-1;
costyup = zeros(T,S,size(ky,2))-1;
costydown = zeros(T,S,size(ky,2))-1;

total = T*S*size(kx,2)*size(ky,2);
for t=1:T
   if(mod(t-1,3)==0)
      disp(sprintf(' Compute norm x cross cost %dx%d/%dx%d\n',t,s,T,S));
   end
   for s=1:S
      %init neightboring mimg
      yup =zeros(H,W,C);
      ydown = zeros(H,W,C);
      xleft = yup;
      xright = zeros(H,W,C);
      % extract neiboring mimg
      if(t~=1)
         yup = squeeze(lf(t-1,s,:,:,:));
      end
      if(t~=T)
         ydown = squeeze(lf(t+1,s,:,:,:));
      end

      if(s~=1)
         xleft = squeeze(lf(t,s-1,:,:,:));
      end
      if(s~=S)
         xright = squeeze(lf(t,s+1,:,:,:));
      end

      % extract subview center
      icenter = squeeze(lf(t,s,centery-hm:centery+hm,centerx-hm:centerx+hm,:));

      for ix=1:size(kx,2)
         costxleft(t,s,ix) = normxcross_rgb(icenter,xleft(centery-hm:centery+hm,kx(ix)-hm:kx(ix)+hm,:));
         costxright(t,s,ix) = normxcross_rgb(icenter,xright(centery-hm:centery+hm,kx(ix)-hm:kx(ix)+hm,:));
      end

      for iy=1:size(ky,2)
         costyup(t,s,iy) = normxcross_rgb(icenter,yup(ky(iy)-hm:ky(iy)+hm,centerx-hm:centerx+hm,:));
         costydown(t,s,iy) = normxcross_rgb(icenter,ydown(ky(iy)-hm:ky(iy)+hm,centerx-hm:centerx+hm,:));
      end
   end
end

disp(sprintf('Finish calculate the norm x cross cost \n'));


disp(sprintf('Merge norm x cross cost \n'));
[dumm cleft] = max(permute(costxleft,[3 1 2]));
cleft=squeeze(cleft);
Cleft=kx(cleft)-centerx;
Cleft = abs(Cleft);

[dumm cright] = max(permute(costxright,[3 1 2]));
cright=squeeze(cright);
Cright=kx(cright)-centerx;
Cright=abs(Cright);
Cleft(:,1) = Cright(:,1);  %cost at border left and right should be corrected.
Cright(:,size(Cright,2)) = Cleft(:,size(Cright,2));

Cx=1/2*(Cleft+Cright);

[dumm cup] = max(permute(costyup,[3 1 2]));
cup=squeeze(cup);
Cup=ky(cup)-centery;
Cup=abs(Cup);

[dumm cdown] = max(permute(costydown,[3 1 2]));
cdown=squeeze(cdown);
Cdown=ky(cdown)-centery;
Cdown=abs(Cdown);
Cup(1,:)=Cdown(1,:); %cost at border up and down should be corrected.
Cdown(size(Cdown,1),:)=Cup(size(Cdown,1),:);

Cy=1/2*(Cup+Cdown);

DepthCost = 1/2*(Cx+Cy);
%after this step the DepthCost is non integer such as k*1/4
%DepthCost is also the diameter of patch.
