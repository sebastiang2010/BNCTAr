%xmp = patient.xstep/2;
%ymp = patient.ystep/2;
%zmp = patient.zstep/2;

%x = patient.xmin+xmp:patient.xstep:patient.xmax-xmp;
%y = patient.ymin+ymp:patient.ystep:patient.ymax-ymp;
%z = patient.zmin+zmp:patient.zstep:patient.zmax-zmp;

% zz=z(10:20);
% yy=y(80:90);
% xx=x(70:80);

%xc=patient.D1c.xc; 
%yc=patient.D1c.yc;
%zc=patient.D1c.zc;


clc 
clear 

step=0.3; 

% imagen 
x=-38.25:step:38.25;
y=-38.25:step:38.25;
z=-19.05:step:19.05;

boro_c=rand(106,91,104);

% mesh 
xc=-15.75:step:15.75;
yc=-13.5:step:13.5;
zc=-15.45:step:15.45;

[X,Y,Z] = ndgrid(x,y,z);
[Xc,Yc,Zc] = ndgrid(xc,yc,zc);


boro_int=interpn(Xc,Yc,Zc,boro_c,X,Y,Z);

%[XX,YY,ZZ] = meshgrid(x,y,zz);
%[XXc,YYc,ZZc] = meshgrid(yc,xc,zc); 

%boro2 = interp3(XXc,YYc,ZZc,boro_c,XX,YY,ZZ,'linear');

step=0.3;

x=-38.25:step:38.25;
y=-38.25:step:38.25;
z=-19.05:step:19.05;

z=z(10:20);
y=y(80:90);
x=x(70:80);



xc=-15.75:step:15.75; 
yc=-13.5:step:13.5;
zc=-15.45:step:15.45;

xc=xc(1:40);
yc=yc(1:40);
zc=zc(1:30);


%boro_c=rand(,,);

% Verificar si step es igual en x y z

    [X,Y,Z] = ndgrid(x,y,z);
    [Xc,Yc,Zc] = ndgrid(xc,yc,zc);

    boro_c=rand(size(X));

    boro_int1=interpn(Xc,Yc,Zc,boro_c,X,Y,Z);

% Reemplazar los valores donde corresponda
for i=1:numel(xc)
    for j=1:numel(yc)
        for k=1:numel(zc)
            if boro_int(i,j,k) > boro_c(i,j,k)
                boro_int(i,j,k) = boro_c(i,j,k);
            end
        end
    end
end

figure(1)
pcolor(boro_int(:,:,15))
figure(2)
pcolor(boro_int1(:,:,15))


% Verificar que la asignación se realizó correctamente
%assert(all(boro_int(:) <= boro_c(:)))

