% Carga los datos de la variable que deseas visualizar
clc 
close all 
I=patient.stack; 

a=unique(I); 

nroi=9; 

img_rois = patient.roi_level(nroi);


% Extrae los datos de la isosuperficie en ind=15
%valor_iso = 93;
%umbral_iso = 0.5;
A=I==img_rois;
A=double(A); 
% for i=1:88
%  A1(:,:,i)=imbinarize(A(:,:,i)); 
% end 
A1=zeros(size(I));

for i=1:88 
borde=bwboundaries(A(:,:,50)); 
A1(borde{1,1},50)=1; 
end 


imshow(A1(:,:,50),[])
    figure
    A=I==106; 
    A=double(A);
    sigma=0.3;
for i=1:88
    imshow(A1(:,:,i),[])
    pause(0.5)

end 

    A1 = imgaussfilt3(A, sigma);
    img_rois = patient.roi_level(i);
    isovalores = isosurface(A1,1);
    isoparche = patch(isovalores);
    set(isoparche,'FaceColor',[0 1 0],'EdgeColor','none');
    daspect([1 1 1]);
    view(3);
    axis tight;
    camlight;
    lighting gouraud;
%end