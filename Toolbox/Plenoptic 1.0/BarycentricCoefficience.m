%Barycentric coefficience 
function f=BarycentricCoefficience(x1,y1,...
                                   x2,y2,...
                                   x3,y3,...
                                   x,y)
A=double([x1 x2 x3; y1 y2 y3; 1 1 1]);
A=inv(A);
f = A*double(transpose([x y 1]));
end