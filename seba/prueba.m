clc 
clear
close all
load('patient.mat')
%R = imref2d([2 2],4,2); 
% revisar como se grafica la dosis 

% grayImage = imread('cameraman.tif');
% imshow(flipud(grayImage), 'XData', [0, 3], 'YData', [0, 3]);
% set(gca, 'YDir', 'normal');
% axis on;
% axis image

Phantom=patient.stack; 







h=plot(rand(1,5),'ButtonDownFcn',@lineCallback);








figure(500)
himg=image(Phantom(:,:,5)); 

% b
n1=patient.beam.n1; 
D1=n1.Matrix; 
D1boro=D1(:,:,:,1); 

n2=patient.beam.n2; 
D2=n2.Matrix; 
D2boro=D2(:,:,:,1); 

Dtotal=patient.Dt.boro; 

figure(1)
hinsow=imshow(D1boro(:,:,44),[]); 
colormap(jet)

figure(2)
imshow(D2boro(:,:,44),[])
colormap(jet)

figure(5)
imshow(Dtotal(:,:,44),[])
colormap(jet)

%dif=D2boro-D1boro; 
%a=sum(dif(:));

figure(101)
imshow(Phantom(:,:,44),[])
axis on 
lim=axis;  

figure(105)
imshow(Phantom(:,:,44),[])
hold on 
h=imshow(Dtotal(:,:,44),[]); 
colormap(jet)
set(h,'AlphaData',0.4)
axis on 

R = imref2d(size(Phantom),0.3,0.3);

xc=6.45:0.3:37.95; 
yc=-26.5:0.3:0.75; 
zc=77.85:0.3:7.05;

function lineCallback(src,~)
   src.Color = rand(1,3);
end


