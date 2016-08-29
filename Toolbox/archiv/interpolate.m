
interpolated_image=[];
[NumberOfRow, NumberOfColumn]=size(view5);
scale=5;
for row=2*scale:NumberOfRow*scale-2*scale
    for column=2*scale:NumberOfColumn*scale-2*scale
        r=row/scale-round(row/scale);
        c=column/scale-round(column/scale);
        sigma=2;
        a1=view1(round(row/scale)-1,round(column/scale)-1);
        a2=view2(round(row/scale)  ,round(column/scale)-1);
        a3=view3(round(row/scale)+1,round(column/scale)-1);
        a4=view4(round(row/scale)-1,round(column/scale)  );
        a5=view5(round(row/scale)  ,round(column/scale)  );
        a6=view6(round(row/scale)+1,round(column/scale)  );
        a7=view6(round(row/scale)-1,round(column/scale)+1);
        a8=view8(round(row/scale)  ,round(column/scale)+1);
        a9=view9(round(row/scale)+1,round(column/scale)+1);
        value_view=[a1 a2 a3 a4 a5 a6 a7 a8 a9];
        interpolated_image(row,column)=sum(GaussianWeigh(r,c,sigma)*transpose(value_view));
    end                                   
end 
%interpolated_image=interpolated_image*20;
imshow(interpolated_image,[])