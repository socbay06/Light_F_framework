function f=disparityMap(viewLeft,viewRight)
disparityRange = [-6 10];
disparityMap = f= disparity(viewCenter,viewLeft,'BlockSize',5,'DisparityRange',disparityRange);
a=imshow(disparityMap, disparityRange);
%colorbar
%imwrite(a,'disparityMap-R.png');

%imshow(imresize(a, 3.0));
end

