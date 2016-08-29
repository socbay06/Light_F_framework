%calculate cubic convolution interpolation kernel
%s= abs distance from point to point ([0,2])
function f=CubicConvolutionKernel(s)
    s=abs(s);
    if s<1
        f= 1.5*s^3- 2.5*s^2+ 1;
    elseif s<2
        f= -.5*s^3+ 2.5*s^2- 4*s +2;
    else
        f= 0;   
end