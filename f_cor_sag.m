function f_cor_sag(patient)

% probar usar imfuse


I=patient.stack;
%s=size(I);

px=patient.px;
py=patient.py;
%pz=patient.nimg; 

switch patient.show
    case 0
        D = [];
    case 1
        D = patient.Df;
    case 2
        D = patient.Dw;
    case 3
        D = patient.Diso;
    case 4
        D = patient.Dt.E(:,:,:,patient.nE); % es el error
end

if patient.show>=1
    %dose_axial=D(:,:,pz);
    dose_sag=squeeze(D(px,:,:));
    %dose_sag=dose_sag';
    dose_cor=squeeze(D(:,py,:));
    %dose_cor=dose_cor';
    %dose_ax=squeeze(D(px,:,:));
    %dose_cor=dose_cor';
end


% xmp = patient.xstep/2;
% ymp = patient.ystep/2;
% zmp = patient.zstep/2;
% 
% x = patient.xmin+xmp:patient.xstep:patient.xmax-xmp;
% y = patient.ymin+ymp:patient.ystep:patient.ymax-ymp;
% z = patient.zmin+zmp:patient.zstep:patient.zmax-zmp;
% 
% xW=[x(1) x(end)]; 
% %yW=[y(1) y(end)]; 
% zW=[z(1) z(end)];

%% axial 
% Iax=I(:,:,pz); 
% 
% R=imref2d(size(Iax));
% R.XWorldLimits=xW; 
% R.YWorldLimits=yW;

%ok=sizesMatch(R,Iax); 
% figure(1)
% cla
% imshow(Iax,R)
% if patient.show>=1
%     hold on
%     % Generar un corte coronal
%     h=imshow(dose_axial,[]);
%     set(h,'AlphaData',0.4)
%     %surf(dose_cor,'FaceAlpha',0.4);
%     %set(patient.dose_surf,'ButtonDownFcn',@imgClick)
%     shading interp;
%     colormap jet
% end
% txt=['Coronal :  ',num2str(pz)];
% title(txt)
% axis square

%% sagittal 
Isag=squeeze(I(px,:,:));
Isag=cat(3,Isag,Isag,Isag);

% Generar un corte sagital
% R1=imref2d(size(Isag));
% R1.XWorldLimits=xW; 
% R1.YWorldLimits=zW;
%ok=sizesMatch(R1,Isag); 

figure(2)
ax=gca;
cla(ax)
image(Isag);
hold on
if patient.show>=1
    % Generar un corte coronal
    %h=imshow(dose_cor,[]);
    %set(h,'AlphaData',0.4)
    surf(dose_sag,'FaceAlpha',0.4);
    %set(patient.dose_surf,'ButtonDownFcn',@imgClick)
    shading interp;
    colormap jet
end
txt=['Sagittal :  ',num2str(px)];
title(txt)
axis square



%% Coronal 
Icor=squeeze(I(:,py,:));
Icor=cat(3,Icor,Icor,Icor); 

% R2=imref2d(size(Icor));
% R2.XWorldLimits=zW; 
% R2.YWorldLimits=xW;
%ok=sizesMatch(R2,Icor); 

figure(3)
ax=gca;
cla(ax)
image(Icor);
hold on
if patient.show>=1
    % Generar un corte coronal
    %h=imshow(dose_cor,[]);
    %set(h,'AlphaData',0.4)
    surf(dose_cor,'FaceAlpha',0.4);
    %set(patient.dose_surf,'ButtonDownFcn',@imgClick)
    shading interp;
    colormap jet
end
txt=['Coronal :  ',num2str(py)];
title(txt)
axis square


end % function





