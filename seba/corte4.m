function f_cor_sag(patient)

% Cargar una imagen DICOM
%info = dicominfo('image.dcm');
close all 
 

I=patient.stack; 
I1=I(:,:,patient.nimg); 
I1 = cat(3,I1,I1,I1);

px=patient.px; 
py=patient.py; 

if patient.show==2; D=patient.Dw; 
end 

Isag=squeeze(I(px,:,:));
Isag=cat(3,Isag,Isag,Isag);

Icor=squeeze(I(:,py,:));
Icor=cat(3,Icor,Icor,Icor);

dose_axial=D(:,:,patient.nimg); 
dose_sag=squeeze(D(px,:,:)); 
dose_cor=squeeze(D(:,py,:)); 


figure(1)
image(I1)
%if(~isfield(patient,'dose_surf'))
hold on 
h=surf(dose_axial,'FaceAlpha',0.4);
%set(patient.dose_surf,'ButtonDownFcn',@imgClick)
shading interp;
colormap jet

% Generar un corte sagital
figure(2)
image(Isag);
hold on 
% Generar un corte coronal
h=surf(dose_sag,'FaceAlpha',0.4);
%set(patient.dose_surf,'ButtonDownFcn',@imgClick)
shading interp;
colormap jet

figure(3)
image(Icor);
hold on 
% Generar un corte coronal
h=surf(dose_cor,'FaceAlpha',0.4);
%set(patient.dose_surf,'ButtonDownFcn',@imgClick)
shading interp;
colormap jet

end 





