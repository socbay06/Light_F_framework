%Gaussian Kernel
%kernel = fspecial('gaussian', 5, 2); surf(kernel);
function w= GaussianWeigh(r,c,sigma)
    %r=.3;
    %c=.3;
    %sigma=1;
    w=[];
    distance_x=[-r-1,-r-1,-r-1;...
                -r  ,-r  ,-r  ;...
                -r+1,-r+1,-r+1];
    distance_y=[-c-1,-c  ,-c+1;...
                -c-1,-c  ,-c+1;...
                -c-1,-c  ,-c+1];
    for i=1:9
        if or(distance_x(i).^2>1,distance_y(i).^2>1)
            w(i)=0;
        else
            w(i)= 1/(2*pi()*sigma^2)*exp(-(distance_x(i).^2+distance_y(i).^2)/(2*sigma^2));
        end
    end
end
