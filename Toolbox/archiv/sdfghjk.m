interpolated_image=[];
[xxx,yyy, NumberOfRow, NumberOfColumn, color]=size(LF);
scale=5;
sigma=.5;
half_windows_size=1;

for row=2*scale:NumberOfRow*scale-2*scale
    for column=2*scale:NumberOfColumn*scale-2*scale
        interpolated_image(row,column)=0;
        r=round(row/scale);
        c=round(column/scale);
        interpolated_pixel_r= (row/scale-round(row/scale))+5;
        interpolated_pixel_c= (column/scale-round(column/scale))+5;
        
        for pixel_r=5-half_windows_size:5+half_windows_size
            for pixel_c=5-half_windows_size:5+half_windows_size
                if or ( abs(pixel_r-interpolated_pixel_r)>half_windows_size,...
                        abs(pixel_c-interpolated_pixel_c)>half_windows_size)
                    continue
                end
                interpolated_image(row,column)=...
                    interpolated_image(row,column)+...
                    1/(2*pi()*sigma^2)*exp(-((pixel_r-interpolated_pixel_r).^2+(pixel_c-interpolated_pixel_c).^2)/(2*sigma^2))*...
                    LF(pixel_c,pixel_r, r+(5-pixel_r), c+(5-pixel_c), 1);
            end
        end
    end
end
imshow(interpolated_image,[])