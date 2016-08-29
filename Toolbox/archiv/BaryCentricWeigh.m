%Barycentric coefficience 

%   1   4   7
%   
%   2   5   8
%
%   3   6   9

function w= BaryCentricWeigh(r,c,Scale)
f=[]; 
w=[];
A=[ round(c/Scale),round(c/Scale), round(c/Scale);...
    round(r/Scale),round(r/Scale), round(r/Scale);... 
    1              1               1];


if  and((round(c/Scale)-c/Scale)>=0,...
        (round(r/Scale)-r/Scale)>0)
    A(2)=A(2)-1;
    A(7)=A(7)+1;
    %4 5 8
    A=inv(A);
    f= A*transpose([c/Scale r/Scale 1]);
    w=[0 0 0 f(1) f(2) 0 0 f(3) 0];
elseif and((round(c/Scale)-c/Scale)>0,...
        (round(r/Scale)-r/Scale)<=0)
    A(5)=A(5)+1;
    A(7)=A(7)+1;
    %5 6 8
    A=inv(A);
    f = A*transpose([c/Scale r/Scale 1]);
    w=[ 0 0 0 0 f(1) f(2) 0 f(3) 0];
elseif and((round(c/Scale)-c/Scale)<=0,...
           (round(r/Scale)-r/Scale)<0)
    A(1)=A(1)-1;
    A(8)=A(8)+1;
    %2 5 6
    A=inv(A);
    f = A*transpose([c/Scale r/Scale 1]);
    w=[ 0 f(1) 0 0 f(2) f(3) 0 0 0];
elseif and((round(c/Scale)-c/Scale)<0,...
           (round(r/Scale)-r/Scale)>=0)
    A(1)=A(1)-1;
    A(5)=A(5)-1;
    %2 4 5
    A=inv(A);
    f = A*transpose([c/Scale r/Scale 1]);
    w=[ 0 f(1) 0 f(2) f(3) 0 0 0 0];
else and((round(c/Scale)-c/Scale)==0,...
        (round(r/Scale)-r/Scale)==0)
    %5
    A=inv(A);
    f = A*transpose([c/Scale r/Scale 1]);
    %w=[ 0 0 0 0 f(1) 0 0 0 0];
    w=[ 0 .25 0 .25 0 .25 0 .25 0];
end

end