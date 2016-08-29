%Barycentric interpolation
% interpolating_point(X) is calculated base on 3 neighbor_pixels, which
% form a matrix.
%
% (x1,y1)        X  
%
%                 o
% (x3,y3)  X            X  (x3,y3)
%
% Value of "x" point: g(x,y)=lambda_1*f(x1,y1)+lambda_2*f(x2,y2)+lambda_3*f(x3,y3);
% where lambda_1,lambda_2,lambda_3 should satisfied system of equations: 
%        x = x1*lambda_1+ x2*lambda_2+ x3*lambda_3
%        y = y1*lambda_1+ y2*lambda_2+ y3*lambda_3
%        1 =    lambda_1+    lambda_2+    lambda_3
Scale=3;
[Smallrow Smallcolumn ]= size(view5);
Bigrow=Scale*Smallrow;
Bigcolumn=Scale*Smallcolumn;
Bigimage=[];
%for co=2*Scale:Bigcolor-2*Scale
    for r=2*Scale:Bigrow-12
        for c=2*Scale:Bigcolumn-12
        a1=view1(round(r/Scale)-1,round(c/Scale)-1);
        a2=view2(round(r/Scale)  ,round(c/Scale)-1);
        a3=view3(round(r/Scale)+1,round(c/Scale)-1);
        a4=view4(round(r/Scale)-1,round(c/Scale)  );
        a5=view5(round(r/Scale)  ,round(c/Scale)  );
        a6=view6(round(r/Scale)+1,round(c/Scale)  );
        a7=view6(round(r/Scale)-1,round(c/Scale)+1);
        a8=view8(round(r/Scale)  ,round(c/Scale)+1);
        a9=view9(round(r/Scale)+1,round(c/Scale)+1);
        value_view=[a1 a2 a3 a4 a5 a6 a7 a8 a9]; 
                                              
            Bigimage(r,c)= sum(BaryCentricWeigh(r,c,Scale)*transpose(value_view));
        end
    end
%end
imshow(Bigimage,[]);
