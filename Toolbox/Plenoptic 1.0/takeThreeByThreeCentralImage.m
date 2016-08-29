% LF matrix DIM = [9 9 381 382 4]
% each microlens image contain 9x9 pixels <=> matrix 9x9
% extract column+row  5  
function extractSubImage(LF)
[NumberofRowMicroImage NumberofColumnMicroImage NumberOfRow NumberOfColumn colorSystem]=size(LF);
    extract_r =[];
    for row=1:NumberOfRow
        temp=[];
        for col=1:NumberOfColumn
            unit=LF([6],[5],row,col,1);
            temp =cat(2,temp, unit);
        end
        extract_r =cat(1,extract_r,temp);
    end
    extract_gr =[]; 
    for row=1:NumberOfRow
        temp=[];
        for col=1:NumberOfColumn
            unit=LF([6],[5],row,col,2);
            temp =cat(2,temp, unit);
        end
        extract_gr =cat(1,extract_gr,temp);
    end
    extract_gb =[];
    for row=1:NumberOfRow
        temp=[];
        for col=1:NumberOfColumn
            unit=LF([6],[5],row,col,3);
            temp =cat(2,temp, unit);
        end
        extract_gb =cat(1,extract_gb,temp);
    end
    extract_b =[];
    for row=1:NumberOfRow
        temp=[];
        for col=1:NumberOfColumn
            unit=LF([6],[5],row,col,4);
            temp =cat(2,temp, unit);
        end
        extract_b =cat(1,extract_b,temp);
    end
    %color analog gain:
    % r     4.375
    % gr    3.2188
    % b     3.7188
    blue=extract_b/1.17;
    green=extract_gr/1.35;
    
    view2=cat(3,extract_r,green,blue);               %output
    %viewFulRGB = demosaic(viewFull_left,'bggr');
    image(view2)% export Image for RGB standard
    imwrite(view2,'view2.png');
end