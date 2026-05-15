%clear 
clc 
close all
%clear 

I=patient.stack; 
A=patient.D1c; 
xstep=0.3; 
ystep=0.3; 
zstep=0.3; 

xmp = patient.xstep/2;
ymp = patient.ystep/2;
zmp = patient.zstep/2;

x = patient.xmin+xmp:patient.xstep:patient.xmax-xmp;
y = patient.ymin+ymp:patient.ystep:patient.ymax-ymp;
z = patient.zmin+zmp:patient.zstep:patient.zmax-zmp;

% xstep_mesh=0.3; 
% ystep_mesh=0.3; 
% zstep_mesh=0.3; 

%xc=-15.75:xstep_mesh:15.75; 
%yc=-13.5:ystep_mesh:13.5;
%zc=-15.45:zstep_mesh:15.45;

xc=A.xc;
yc=A.yc; 
zc=A.zc; 

[Xc,Yc,Zc] = ndgrid(xc,yc,zc);
[X,Y,Z] = ndgrid(x,y,z);

boro_c=rand(size(Xc));

s1=size(Xc); 
s = size(X);

boro_int= nan(s);

x_idx = find(xc(1)== x) ;% & xc== x(end));
y_idx = find(yc(1)== y);% & yc <= y(end));
z_idx = find(zc(1)== z); %& zc <= z(end));

%boro_int(x_idx:x_idx+s1(1),y_idx:y_idx+s1(2) ,z_idx:z_idx+s1(3)) = boro_c;





