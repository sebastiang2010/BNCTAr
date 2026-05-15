function D1c=f_lectura_phits_neutron 

clc
%clear
%close all

format longG

%% agregar la lectura de nx ny nz
%% si part es neutron agrgar que hay que leer 3 mset

%filesave='KERMA_neutrones.mat'; 

currentdirectory=pwd;
%if isempty(file)
tipoarchivo='*.out';
[archivo,directorio]=uigetfile(tipoarchivo,'Select File');
[~,msg]=fopen(fullfile(directorio,archivo));
if ~isempty(msg)
    disp(' ')
    disp(msg)
    return
end
file=fullfile(directorio,archivo);

%tic
lines=readlines(file);
%time2=toc;



txt=' newpage:';
txt1=' part ';
txt2=' nx =';
txt3=' ny =';
txt4=' nz =';
%txt5= ' mset =';
txt10{1,1}=' xmin ='; 
txt10{1,2}=' xmax = ';  
txt10{1,3}=' xdel ='; 
txt10{1,4}=' ymin ='; 
txt10{1,5}=' ymax =';  
txt10{1,6}=' ydel =';
txt10{1,7}=' zmin ='; 
txt10{1,8}=' zmax =';  
txt10{1,9}=' zdel ='; 


for i = 1:length(txt10)
    
    b = find(contains(lines,txt10{1,i}));   
    p = strtrim(lines(b(1)));
    p = strsplit(p, '=');
    
    valor_str = strtrim(p{2});
    
    num = regexp(valor_str, '[-+]?\d*\.?\d+', 'match');
    datos_mesh{1,i} = str2double(num{1});
    
end
xmin = datos_mesh{1,1}; 
xmax = datos_mesh{1,2};
ymin = datos_mesh{1,4};
ymax = datos_mesh{1,5}; 
zmin = datos_mesh{1,7} ; 
zmax = datos_mesh{1,8}; 
 

xstep = datos_mesh{1,3};
ystep = datos_mesh{1,6}; 
zstep = datos_mesh{1,9}; 

xmp = xstep/2;
ymp = ystep/2;
zmp = zstep/2;

xc = xmin+xmp:xstep:xmax-xmp;
yc = ymin+ymp:ystep:ymax-ymp;
zc = zmin+zmp:zstep:zmax-zmp;



%%
% busco particula
b=find(contains(lines,txt1));
p = strtrim(lines(b(1)));
p = strsplit(p, '=');
part=strtrim(p{2});
clear p b

if strcmpi(part, 'neutron')
    checkmset='true';
    if checkmset==false
        disp(' NO ES NEUTRON')
         
    end    
end


% busco nx, ny, nz
c=find(contains(lines,txt2));
p = strtrim(lines(c(1)));
nx = regexp(p, '\d+', 'match');
nx=str2double(nx);

c=find(contains(lines,txt3));
p = strtrim(lines(c(1)));
ny = regexp(p, '\d+', 'match');
ny=str2double(ny);

c=find(contains(lines,txt4));
p = strtrim(lines(c(1)));
nz = regexp(p, '\d+', 'match');
nz=str2double(nz);

clear c p

% busco los newpage
a=find(contains(lines,txt));
a1=find(lines=='#newpage:'); % la primera es diferente


a=[a1,a']; 

clear a1 
m=round(nx*ny/10);

%nz=nz-1;
salto=25; %desde newpage a donde empieza la matriz
m1=m-1; % tamaño de la lectura
saltomset=4;
% veo que mset es

C=ones(nx,ny,nz,3); 

for g=1:3


    %p=lines(a1+saltomset);
    %p = strtrim(p);
    %p = strsplit(p, '=');
    %nmset=strtrim(p{2});
    %nmset=str2double(nmset);
    
    if g==1
        % esto es diferente porque el primero tiene #newpage
        nmset=g;
        n=1; 
        for k=1:nz

            %p=lines(a(k)+saltomset);
            %p = strtrim(p);
            %p = strsplit(p, '=');
            %nmset=strtrim(p{2});
            %nmset=str2double(nmset);

            A=lines(a(k)+salto:a(k)+salto+m1);


            B10 = cell2mat(cellfun(@(x) str2double(strsplit(strtrim(x))), A(1:end-1), 'UniformOutput', false));
            B11 =  cell2mat(cellfun(@(x) str2double(strsplit(strtrim(x))), A(end,:), 'UniformOutput', false));
            B12 = NaN(1,10);
            B12(1:length(B11))=B11;

            B13=[B10;B12];

            if mod(nx*ny,10)==0
                B=reshape(B13',nx,ny);
                B1=B';
            else
                A1=B13';
                A2=A1(:);
                %B=zeros(nx,ny);
                %nmset=lines(a(k)+saltomset);
                for j=1:ny
                    o=A2(1:nx);
                    B(:,j)=o';
                    A2(1:nx)=[];
                end
                B1=B;
            end

            C(:,:,n,nmset)=B1;
            n=n+1; 
        end
    end
    if g==2
        nmset=2;
        n=1;
        for k=nz+1:2*nz

            %p=lines(a(k)+saltomset);
            %p = strtrim(p);
            %p = strsplit(p, '=');
            %nmset=strtrim(p{2});
            %nmset=str2double(nmset);

            A=lines(a(k)+salto:a(k)+salto+m1);


            B10 = cell2mat(cellfun(@(x) str2double(strsplit(strtrim(x))), A(1:end-1), 'UniformOutput', false));
            B11 =  cell2mat(cellfun(@(x) str2double(strsplit(strtrim(x))), A(end,:), 'UniformOutput', false));
            B12 = NaN(1,10);
            B12(1:length(B11))=B11;

            B13=[B10;B12];

            if mod(nx*ny,10)==0
                B=reshape(B13',nx,ny);
                B1=B';
            else
                A1=B13';
                A2=A1(:);
                %B=zeros(nx,ny);
                %nmset=lines(a(k)+saltomset);
                for j=1:ny
                    o=A2(1:nx);
                    B(:,j)=o';
                    A2(1:nx)=[];
                end
                B1=B;
            end

            C(:,:,n,nmset)=B1;
            n=n+1;
        end
    end

    if g==3
        nmset=3;
        n=1;
        
        for k=2*nz+1:3*nz

            %p=lines(a(k)+saltomset);
            %p = strtrim(p);
            %p = strsplit(p, '=');
            %nmset=strtrim(p{2});
            %nmset=str2double(nmset);

            A=lines(a(k)+salto:a(k)+salto+m1);


            B10 = cell2mat(cellfun(@(x) str2double(strsplit(strtrim(x))), A(1:end-1), 'UniformOutput', false));
            B11 =  cell2mat(cellfun(@(x) str2double(strsplit(strtrim(x))), A(end,:), 'UniformOutput', false));
            B12 = NaN(1,10);
            B12(1:length(B11))=B11;

            B13=[B10;B12];

            if mod(nx*ny,10)==0
                B=reshape(B13',nx,ny);
                B1=B';
            else
                A1=B13';
                A2=A1(:);
                B=zeros(nx,ny);
                %nmset=lines(a(k)+saltomset);
                for j=1:ny
                    o=A2(1:nx);
                    B(:,j)=o';
                    A2(1:nx)=[];
                end
                B1=B;
            end
             
            if k==60;disp(B1);end
            

            C(:,:,n,nmset)=B1;
            n=n+1;
        end
    end
     
    
end 

%%
%for i=1:3
%    B = flip(C(:,:,:, i)); 
%    C1(:,:,:,i)=B;
%end 
%C=[]; 
%C=C1; 
%clear C1


%%
% figure(105)
% for i=1:nz
%     imshow(C(:,:,i,1),[])
%     colormap(jet)
%     pause(0.1)
% end
% 
% figure(106)
% for i=1:nz
%     imshow(C(:,:,i,2),[])
%     colormap(jet)
%     pause(0.1)
% end
% 
% figure(107)
% for i=1:nz
%     imshow(C(:,:,i,3),[])
%     colormap(jet)
%     pause(0.1)
% end
% 
% 
% figure(800)
% imshow(C(:,:,59,1),[])
% colormap(jet)

%Dn=C; 


%% 
% Dosis.neutron=Dn; 
% Dosis.size=[nx ny nz]; 
% Dosis.voxel=[0.3 0.3 0.3]; 
% Dosis.x=[-15.9 15.9]; 
% Dosis.y=[-13.65 13.65]; 
% Dosis.z=[-15.9 15.9]; 

for i=1:3
    A=C(:,:,:,i);
    A1=A; 
    %A1=permute(A, [3 1 2]);
    if i==1 
        D1c.boro=A1;
    end 
    if i==2
        D1c.fast=A1; 
    end 
    if i==3      
        D1c.thn=A1;
    end 
 end 



D1c.xc=xc; 
D1c.yc=yc; 
D1c.zc=zc; 

cd(currentdirectory);

