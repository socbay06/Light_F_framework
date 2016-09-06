% 4gb channel norm cross correlation
% input: 2 same size matrix.
function out=normxcross_rgb(mat1,mat2)
   out = normxcross_single(mat1(:,:,1),mat2(:,:,1));
   out = out + normxcross_single(mat1(:,:,2),mat2(:,:,2));
   out = out + normxcross_single(mat1(:,:,3),mat2(:,:,3));
   out = 1/3*out;

end
