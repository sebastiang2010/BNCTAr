k%clear 
clc 
%P=patient.stack; 

% %%  funcion handle 
% f=@sin;  
% f(pi/4) 
% sin(pi/4)
% % 
% g=@(x) x^2; 
% r=@(x,y) sqrt(x.^2+y^2); 
% 
% g(3)
% r(3,4)
% 
% % analizar 
% % timeit
% 
% 
% %% 
% nivel_gris=unique(P);
% IND_roi = P == reshape(nivel_gris, [1 1 1 numel(nivel_gris)]);
% 
% %% meshgrid 
% [X,Y]=meshgrid(1:3,10:14);

%% 
xmp = patient.xstep/2;
ymp = patient.ystep/2;
zmp = patient.zstep/2;

x = patient.xmin+xmp:patient.xstep:patient.xmax-xmp;
y = patient.ymin+ymp:patient.ystep:patient.ymax-ymp;
z = patient.zmin+zmp:patient.zstep:patient.zmax-zmp;

[X,Y,Z] = ndgrid(x,y,z);


% Visualizar la grilla en 3D
% figure;
% scatter3(Y(:), Z(:)) %, 10, 'filled');
% title('Grilla en 3D generada con ngrid');
% xlabel('X'); ylabel('Y'); zlabel('Z');

figure(100)
scatter(Y(:), X(:),0.9,'MarkerFaceAlpha', 0.5);
title('Grilla en 2D generada con ngrid');
xlabel('X'); ylabel('Y');
hold on 
I=patient.stack; 
h=imshow(I(:,:,53),[]);
set(h, 'AlphaData', 0.5);


[Xc,Yc,Zc] = ndgrid(patient.D1c.xc,patient.D1c.yc,patient.D1c.zc);
figure(101)
scatter(Yc(:), Xc(:),0.9,'MarkerFaceAlpha', 0.5);
title('Grilla en 2D generada con ngrid');
xlabel('X'); ylabel('Y');
hold on 
I=patient.stack; 
h=imshow(I(:,:,53),[]);
set(h, 'AlphaData', 0.5);