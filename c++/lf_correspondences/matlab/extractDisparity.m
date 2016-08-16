FILE='lab.h5';
DATAROOT='/flow/';


left=h5read(FILE, strcat(DATAROOT,'left'));
right=h5read(FILE, strcat(DATAROOT,'right'));
up=h5read(FILE, strcat(DATAROOT,'up'));
down=h5read(FILE, strcat(DATAROOT,'down'));

HEIGHT=size(left,1);
WIDTH=size(left,2);

limg=zeros(HEIGHT,WIDTH,3);
rimg=zeros(HEIGHT,WIDTH,3);
uimg=zeros(HEIGHT,WIDTH,3);
dimg=zeros(HEIGHT,WIDTH,3);

tmpp = left(:,:,1);
tmpp(tmpp<0)=0;
limg(:,:,1)=abs(tmpp)/max(tmpp(:));

tmpm = left(:,:,1);
tmpm(tmpm>0)=0;
limg(:,:,2)=abs(tmpm/min(tmpm(:)));

tmpp = right(:,:,1);
tmpp(tmpp<0)=0;
rimg(:,:,1)=abs(tmpp)/max(tmpp(:));

tmpm = right(:,:,1);
tmpm(tmpm>0)=0;
rimg(:,:,2)=abs(tmpm/min(tmpm(:)));


tmpp = up(:,:,1);
tmpp(tmpp<0)=0;
uimg(:,:,1)=abs(tmpp)/max(tmpp(:));

tmpm = up(:,:,1);
tmpm(tmpm>0)=0;
uimg(:,:,2)=abs(tmpm/min(tmpm(:)));


tmpp = down(:,:,1);
tmpp(tmpp<0)=0;
dimg(:,:,1)=abs(tmpp)/max(tmpp(:));

tmpm = down(:,:,1);
tmpm(tmpm>0)=0;
dimg(:,:,2)=abs(tmpm/min(tmpm(:)));

subplot(2,2,1); imshow(limg);
title('left');
subplot(2,2,2); imshow(rimg);
title('right');
subplot(2,2,3); imshow(uimg);
title('up');
subplot(2,2,4); imshow(dimg);
title('down');
