% implement blending algorithm of Georgiev.
% Georgiev, T., & Lumsdaine, A. (2010). Focused plenoptic camera and rendering. Journal of Electronic Imaging, 19(2), 021106. http://doi.org/10.1117/1.3442712
lf=hdf5read('focused.flower.h5','/lf/focused/flower');

[T,S,H,W,C] = size(lf);

M=15;
hM=floor((M-1)/2);
centerx = floor((W-1)/2)+1;
centery = floor((H-1)/2)+1;

limx = floor((floor(W/M)-1)/2);
limy = floor((floor(H/M)-1)/2);
% if lim=0 no blending at all.

view = zeros(T,S,M,M,3);
view0= lf(:,:,centery-hM:centery+hM, centerx-hM:centerx+hM,:);

blenddata=zeros(limy,limx,T,S,M,M,3);
weight = fspecial('average',2*limx+1);
for j=-limy:limy
   for i=-limx:limx
      shifty = centery+(M)*j;
      shiftx = centerx+(M)*i;
      for t=1:T
         for s=1:S
            dshifty = shifty;
            dshiftx = shiftx;
            dt =t+j;
            ds =s+i;
            if(dt<1 || ds<1 || dt>T || ds>S )
               % when there is no data to blend, set it to the center image.
               dt=t;
               ds=s;
               dshifty = centery;
               dshiftx = centerx;
            end
            blenddata(j+limy+1,i+limx+1,t,s,:,:,:) = weight(j+limy+1,i+limx+1)*lf(dt,ds,dshifty-hM:dshifty+hM,dshiftx-hM:dshiftx+hM,:);
         end
      end
   end
end


% blenddata limy x limx x T x S x H x W x C
BData = sum(blenddata,1);
BData = sum(BData,2);
BData = squeeze(BData)/sum(weight(:));
view = BData;

view=flipdim(view,4);
view=flipdim(view,3);

render = reshape(permute(view,[3 1 4 2 5]),[M*T M*S 3]);
figure; imagesc(render)
