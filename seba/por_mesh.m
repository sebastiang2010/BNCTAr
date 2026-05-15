clc 
%clear 
close all 
%clear

%falta aire porque como NaN=0 lo considera aire
% pensar ponerle otro valor 

I=patient.stack; 
Dt=patient.Dw; 
s=size(Dt); 

IND=isnan(Dt);
smesh=size(IND); 

mesh=zeros(s); 
mesh(~IND)=1;  

%figure(1)
%imshow(mesh(:,:,30),[])

%% calculo de voluem que deja afuera la mesh 
P2=double(I).*mesh; 

roi_name={'PTV_BNCT'
    'GTV_3mm'
    'BONE'
    'LEFT_LUNG'
    'BODY'
    'RT_LUNG'
    'RT_INNER_EAR'
    'RT_PAROTID'
    'RT_EYE'
    'LT_EYE'
    'LT_INNER_EAR'
    'LT_OPTIC_NERVE'
    'GTV70'
    'SPINAL_CORD'
    'BRAIN_STEM'
    'LT_LENS'
    'BRAIN'
    'RTOPTIC_NERVE'
    'RT_LENS'
    'LT_PAROTID'
    'AIRE'};

roi_name=flip(roi_name);

roi=unique(I); 
nroi=length(roi); 

meshinfo=struc([]);

%meshinfo.roi{}

%Ojo el aire esta mal hay que sacar lo de la mesh 
for i=2:nroi
    meshinfo.roi{i,1}=roi_name(i); 

    IND=I==roi(i); %imagen real 
    a=sum(IND(:)); 
    %meshinfo.roi{i,2}=a;
    
    IND1=P2==roi(i); %imagen mesh
    b=sum(IND1(:)); 
    %meshinfo.roi{i,3}=b; 
    
    %outside grid
    meshinfo.roi{i,2}=100-b/a*100; 

% 
%         if i==17
%             for j=1:88
%                 subplot(2,2,1)
%                 imshow(IND(:,:,j),[])
%                 subplot(2,2,2)
%                 imshow(IND1(:,:,j),[])
%                 subplot(2,2,3)
%                 imshow(IND(:,:,j).*IND1(:,:,j))
%               pause
%             end
%         end

    
end           








