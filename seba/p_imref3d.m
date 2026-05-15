clc 


R=imref3d(size(I));
xW=[x(1) x(end)]; 
yW=[y(1) y(end)]; 
zW=[z(1) z(end)];
%TF = contains(R,xWorld,yWorld,zWorld)

R.XWorldLimits=xW; 
R.YWorldLimits=yW;
R.ZWorldLimits=zW;

R1=imref2d(size(I)); 


R.XWorldLimits=xW; 
R.YWorldLimits=yW;
R.ZWorldLimits=zW;
%[xIntrinsic,yIntrinsic,zIntrinsic] = worldToIntrinsic(R,xW,yW,zW);

R1.XWorldLimits=xW; 
R1.YWorldLimits=yW;
%R.ZWorldLimits=zW;


ok = sizesMatch(R,I); 
ok1 = sizesMatch(R1,I); 
%[i, j, k] = worldToSubscript(R,xW,yW,zW)

figure(101)
imshow(I(:,:,70),R1) 
axis on 


xLimits = xlim(gca);