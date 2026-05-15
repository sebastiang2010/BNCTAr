xmp = patient.xstep/2;
ymp = patient.ystep/2;
zmp = patient.zstep/2;

x = patient.xmin+xmp:patient.xstep:patient.xmax-xmp;
y = patient.ymin+ymp:patient.ystep:patient.ymax-ymp;
z = patient.zmin+zmp:patient.zstep:patient.zmax-zmp;

% zz=z(10:20);
% yy=y(80:90);
% xx=x(70:80);

xc=patient.D1c.xc; 
yc=patient.D1c.yc;
zc=patient.D1c.zc;




boro_c=rand(106,91,104);



[X,Y,Z] = ndgrid(x,y,z);
[Xc,Yc,Zc] = ndgrid(xc,yc,zc);


boro = interpn(Xc,Yc,Zc,boro_c,X,Y,Z);

% [XX,YY,ZZ] = meshgrid(x,y,zz);
% [XXc,YYc,ZZc] = meshgrid(yc,xc,zc); 
% 
% boro2 = interp3(XXc,YYc,ZZc,boro_c,XX,YY,ZZ,'linear');


xc=patient.D1c.xc; 
yc=patient.D1c.yc;
zc=patient.D1c.zc;



[Xc,Yc,Zc] = ndgrid(xc,yc,zc,'linear');
