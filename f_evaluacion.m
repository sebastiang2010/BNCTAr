function  component= f_evaluacion(patient)
stack=patient.stack; 

E=patient.Dt.E;
E=E.*100; 

%archivo='C:\MAT\info.htm'; 

name_compoment={'boro','thn','fast','gamma'};

D=zeros(size(E));
D(:,:,:,1)=patient.Dt.boro;
D(:,:,:,2)=patient.Dt.thn;
D(:,:,:,3)=patient.Dt.fast;
D(:,:,:,4)=patient.Dt.g;

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

%[filename,pathname] = uiputfile('*.xls', 'Save info as...');
%archivo=[pathname,filename];
for i=1:4 
    D1=D(:,:,:,i);
    D(:,:,:,i)=D(:,:,:,i)./max(D1(:)).*100; 
end 

for j=1:4
   for i=2:length(niveles)   
    IND_dose(:,:,:,i,j)=(D(:,:,:,j)>niveles(i-1) & D(:,:,:,j)<=i); 
  end     
end 

for i=size(niveles_gris,1)
    IND_roi(:,:,:,i)=stack==nivel_gris(i);
end 



for j=1:4
    E1=E(:,:,:,j);
    %D1=D(:,:,:,j);
    IND_dose1=IND_dose(:,:,:,:,j);
    for i=2:size(nivel_gris,1) % el primero es aire no me interesa
        for k=1:length(niveles)

            IND_roi1=IND_roi(:,:,:,i);
            mask=IND_dose1(:,:,:,:,k).*IND_roi1;


            %IND=stack==nivel_gris(i);
            a=E1(mask);
            max_e=max(a(:));
            min_e=min(a(:));
            desv_e=std(a(:));




            component(j).info{1,i}='Componente'; %componet
            component(j).info{2,i}=name_compoment{1,j}; %componet
            component(j).info{3,i}='Nivel de gris';
            component(j).info{4,i}=nivel_gris(i); %nivel de gris
            component(j).info{5,i}='N ROI';
            component(j).info{6,i}=roi(i);
            component(j).info{7,i}='ROI name';
            component(j).info{8,i}=roi_name{i,1};
            component(j).info{9,i}='Nivel de dosis';
            component(j).info{9,i}=nu2srt(k);
            component(j).info{10,i}=min_e;
            component(j).info{11,i}=max_e;
            component(j).info{12,i}=desv_e;
            %component(j).info{13,i}='------';
            %component(j).info{14,i}=numstr()
            %component(j).info{14,i}=min_d;
            %component(j).info{15,i}=max_d;
            %component(j).info{16,i}=desv_d;

            %writestruct(component,archivo,'FileType','xml')

        end
    end

end