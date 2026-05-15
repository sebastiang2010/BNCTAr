function [component2,component]=f_error(patient)
clc 
%clear 
%close all 

% Poner a zero las variables 

maximoglobalaire=1; %1 el maximo en el aire
                    %0 el maximo en el cuerpo


v=[0 10 20 100]; % rango en el cual analizo el porcentaje de error                    

tvoxel=[0.3 0.3 0.3 ]; %tamaño de la mesh cm3 

stack=patient.stack; 

INDaire=stack==0; 

E=patient.Dt.E;
E=E.*100; 

name_compoment={'boro','thn','fast','gamma'};

%saltos_dosis=5;  
%dose_levels=0:saltos_dosis:100; 
dose_levels=[0 10 100]; % lo que solicita Hanna

%% dosis por componete 
D=zeros(size(E));
D(:,:,:,1)=patient.Dt.boro;
D(:,:,:,2)=patient.Dt.thn;
D(:,:,:,3)=patient.Dt.fast;
D(:,:,:,4)=patient.Dt.g;

s=size(E(:,:,:,1)); 

if maximoglobalaire==0 
   D=D.*(~INDaire); 
end 

E=E.*(~INDaire);

nivel_gris=unique(stack);

roi=1:size(nivel_gris,1);
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

component=struc([]);
component2=struc([]);

a=zeros(1,4);
for i=1:4 
    D1=D(:,:,:,i);
    E1=E(:,:,:,i);
    a(i)=max(E1(:));
    D(:,:,:,i)=D(:,:,:,i)./max(D1(:)).*100; 
end 
IND_nan=isnan(D); 
D(IND_nan)=-20; 
IND_nan=isnan(E); 
E(IND_nan)=-20; 
clear D1

% ind=find(E(:,:,:,1)==a(1)); 
% [x ,y, z]=ind2sub(s,ind);   

tic 
for j=1:4
   for i=1:length(dose_levels)-1  
       IND_dose(:,:,:,i,j)=(D(:,:,:,j)>dose_levels(i) & D(:,:,:,j)<=dose_levels(i+1)); 
  end     
end 
time=toc; 
% igual al de arriba 
% IND_dose = (D > dose_levels(1:end-1)) & (D <= dose_levels(2:end));
% IND_dose = cat(4, false(size(IND_dose(:,:,:,1))), IND_dose);
% IND_dose = permute(IND_dose, [1 2 3 5 4]);

for i=1:size(nivel_gris,1)
     IND_roi(:,:,:,i)=stack==nivel_gris(i);
 end 
%igual al for 
%IND_roi = P == reshape(nivel_gris, [1 1 1 numel(nivel_gris)]);

Phantom=~INDaire; 
%ntotal=sum(Phantom(:)); 

% version Seba
for j=1:4
    E1=E(:,:,:,j);
    IND_dose1=IND_dose(:,:,:,:,j);
    
    for i=2:size(nivel_gris,1) % el primero es aire no me interesa
        
        A=IND_roi(:,:,:,i);

        for k=1:length(dose_levels)-1

    
        

        B=IND_dose1(:,:,:,k); 
        C=logical(A.*B);
        

        b=E1(C);          
            
            %b=[]; 
            a=sum(C(:)); 
            if a>=1
               porcentaje=a/sum(A(:)); 
               porcentaje=porcentaje*100; 
               max_e=max(b(:));
               min_e=min(b(:));
               mean_e=mean(b(:)); 
               desv_e=std(b(:));             
            else 
                porcentaje=0; 
                max_e=[]; 
                min_e=[]; 
                desv_e=[]; 
                mean_e=[];
             end 

           
            component(j).roi(i).info{1,k}='Componente'; %componet
            component(j).roi(i).info{2,k}=name_compoment{1,j}; %componet
            component(j).roi(i).info{3,k}='Nivel de gris';
            component(j).roi(i).info{4,k}=nivel_gris(i); %nivel de gris
            component(j).roi(i).info{5,k}='N ROI';
            component(j).roi(i).info{6,k}=roi(i);
            component(j).roi(i).info{7,k}='ROI name';
            component(j).roi(i).info{8,k}=roi_name{i,1};
            component(j).roi(i).info{9,k}='Dosis (%)';
            txt=[num2str(dose_levels(k)),'-',num2str(dose_levels(k+1))];
            component(j).roi(i).info{10,k}=txt; 
            component(j).roi(i).info{11,k}='Roi (%)';  
            component(j).roi(i).info{12,k}=porcentaje;
            component(j).roi(i).info{13,k}='Min';
            component(j).roi(i).info{14,k}=min_e;
            component(j).roi(i).info{15,k}='Max';
            component(j).roi(i).info{16,k}=max_e;
            component(j).roi(i).info{17,k}='Mean';
            component(j).roi(i).info{18,k}=mean_e;  
            component(j).roi(i).info{19,k}='Desv';
            component(j).roi(i).info{20,k}=desv_e;              
        end %k         
    end %i    
end % j 

%% version sara 
Phantom=~INDaire; 
ntotal=sum(Phantom(:)); 

% figure(100)
% imshow(Phantom(:,:,70),[])

k=1; 
component2.info{1,k}=[];
component2.info{2,k}=[];
component2.info{3,k}=[]; 
component2.info{4,k}=[];
component2.info{5,k}=[]; 
component2.info{6,k}='min';
component2.info{7,k}='max';
component2.info{8,k}='mean'; 
component2.info{9,k}='desv';

for j=1:4
   IND_dose1=IND_dose(:,:,:,:,j); 
   E1=E(:,:,:,j);
   
   component2(j).info{1,1}=[];
   component2(j).info{2,1}=[];
   component2(j).info{3,1}=[];
   component2(j).info{4,1}=[];
   component2(j).info{5,1}=[];
   component2(j).info{6,1}='min';
   component2(j).info{7,1}='max';
   component2(j).info{8,1}='mean';
   component2(j).info{9,1}='desv';
  
   for k=1:length(dose_levels)-1
                       
        B=IND_dose1(:,:,:,k); 
        C=logical(Phantom.*B);

        b=E1(C);

        max_e=max(b(:));
        min_e=min(b(:));
        mean_e=mean(b(:)); 
        desv_e=std(b(:)); 
        
        txt1=[num2str(dose_levels(k)),'-',num2str(dose_levels(k+1))];    
        porcentaje_error=cell(6,size(v,2)-1);% inicializo la celda
        nt=length(b); 
        for i=1:size(v,2)-1
            a=(b>v(i) & b<=v(i+1)); 
            n=sum(a(:)); 
            porcentaje_error{1,i}='Rango de Dosis'; 
            porcentaje_error{2,i}=txt1; 
            porcentaje_error{3,i}='Rango Error (%)';
            txt=[num2str(v(i)),'-',num2str(v(i+1))];
            porcentaje_error{4,i}=txt;
            porcentaje_error{5,i}='Volumen (%)';
            porcentaje_error{6,i}=n/nt*100; 
            porcentaje_error{7,i}='Min';
            porcentaje_error{8,i}=min(b(a));
            porcentaje_error{9,i}='Max';
            porcentaje_error{10,i}=max(b(a));
            porcentaje_error{11,i}='Desv'; 
            porcentaje_error{12,i}=std(b(a)); 
        end 

        
        
        component2(j).info{1,k+1}='Componente'; %componet
        component2(j).info{2,k+1}=name_compoment{1,j}; %componet
        component2(j).info{3,k+1}='Dosis (%)';
        txt=[num2str(dose_levels(k)),'-',num2str(dose_levels(k+1))];
        component2(j).info{4,k+1}=txt;
        component2(j).info{5,k+1}='Error (%)';
        component2(j).info{6,k+1}=min_e;
        component2(j).info{7,k+1}=max_e;
        component2(j).info{8,k+1}=mean_e;
        component2(j).info{9,k+1}=desv_e;
        component2(j).info{10,k+1}='N voxels';
        component2(j).info{11,k+1}=nt; 
        component2(j).info{12,k+1}='Vomuen (cm3)'; 
        component2(j).info{13,k+1}=nt*prod(tvoxel); 
        component2(j).info{14,k+1}='Distr Error'; 
        component2(j).info{15,k+1}=porcentaje_error;

   end 
end 







