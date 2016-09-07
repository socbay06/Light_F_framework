% render the all-in-focused view with Georgiev medthod.
% render final image with best cost.
% input:
% LF = [T,S,H,W,C], T,S: microlens coordinate. H,W: m image resolution.
% PatchSizeMap=[T,S,P];
% view_h,view_w: row,col shift of view.
% target_patch: will scalling to fit this size.
function render=render_aif_georgiev(PatchSizeMap,LF,target_patch, view_h,view_w)

[T,S,H,W,C] = size(LF);

centerx=floor((W-1)/2)+1+view_w;
centery=floor((H-1)/2)+1+view_h;

view0 = zeros(T,S,target_patch,target_patch,3);
for t=1:T
   for s=1:S
      dsize=PatchSizeMap(t,s)*4; % patch size need to be an integer.
      radius = max(floor((dsize-1)/2),1);
      view =squeeze(LF(t,s,:,:,:));
      view4=imresize(view,4);% scale to 4x in oder to keep up with patch size.

      patch = imrotate(view4(4*centery-radius:4*centery+radius,4*centerx-radius:4*centerx+radius,:),180);
      view0(t,s,:,:,:) = imresize(patch,[target_patch target_patch],'bicubic');

   end
end

render = reshape(permute(view0,[ 3 1 4 2 5]),[ T*target_patch S*target_patch 3]);

end
