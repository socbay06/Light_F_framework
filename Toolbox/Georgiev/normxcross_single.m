% single channel norm cross correlation
% input: 2 same size matrix.
function out=normxcross_single(mat1,mat2)
   mat1=squeeze(mat1);
   mat2=squeeze(mat2);
   [H W] = size(mat1);

   mean1 = sum(mat1(:))/numel(mat1);
   vec1=reshape(mat1-mean1,[H*W 1]);
   mean2 = sum(mat2(:))/numel(mat2);
   vec2=reshape(mat2-mean2,[H*W 1]);

   norm1 = norm(vec1);
   norm2 = norm(vec2);

   if(norm1==0 || norm2==0)
      out = -1;
   else
      mulit = vec1.*vec2;
      out = sum(mulit(:))/norm1/norm2;
   end
end
