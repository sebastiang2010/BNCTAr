clc 
clear
close all
 load('patient.mat')
 
 p=patient; 
 clear patient 

%load(fullfile(toolboxdir('images'),'imdata','BrainMRILabeled','images','vol_001.mat'));

D=p.Dt.boro; 
Phantom=p.stack; 

cmap=jet(64); 

figure(200)
imshow(D(:,:,44),[])
colormap jet 


IND=isnan(D); 
D1=D; 
D1(IND)=0; 

figure(102)
s2 = orthosliceViewer(Phantom);% ,'Colormap',cmap); 

figure(103)
sliceViewer(Phantom)


figure(105)
sliceViewer(D1)
colormap jet 

figure(102)
s2 = orthosliceViewer(D1,'Colormap',cmap); 