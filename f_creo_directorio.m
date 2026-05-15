function directorio=f_creo_directorio

if ispc==1
    currentdirectory=pwd;
    addpath(currentdirectory);
    %
    directorio='c:/MAT/BNCTAr/';
    directorio1='d:/MAT/BNCTAr/';
    
    direct=dir(directorio);
    n=size(direct,1);
    
    direct=dir(directorio1);
    n1=size(direct,1);
    %%%%%%%%%%%%%%
    if n>0
        %direct=directorio;
        return
    end
    if n1>0
        directorio=directorio1;
        return
    end
    
   list={'c:/MAT/BNCTAr/','d:/MAT/BNCTAr/'};
   [a,tf] = listdlg('PromptString','Ingrese el directorio','ListString',list,'SelectionMode','single');

    
    if a==1
        
        directorio='c:/MAT/BNCTAr/';
        [s,mess,messid]=mkdir(directorio);
        if s==1
            %clc
            %disp('   ')
            %disp('Se creo el directorio: ')
            %disp(directorio)
            txt=['The directory is created:',directorio]; 
            f = msgbox(txt);
        else
            %disp('   ');
            %disp(mess)
            f = msgbox(mess);
            %dips('   ')
            %dips(messid)
            f = msgbox(messid);
        end
    else
        directorio='d:/MAT/BNCTAr/';
        [s,mess,messid]=mkdir(directorio);
        if s==1
            txt=['The directory is created: ',directorio]; 
            f = msgbox(txt);
        else
            %disp('   ');
            %disp(mess)
            f = msgbox(mess);
            %dips('   ')
            %dips(messid)
            f = msgbox(messid);
        end
    end
    
end

%% caso Unix
if isunix==1
    currentdirectory=pwd;
    addpath(currentdirectory);
    %
    directorio='/Home/MAT/3Dosim';
    
    direct=dir(directorio);
    n=size(direct,1);
    
    %%%%%%%%%%%%%%
    if n>0
        %direct=directorio;
        return
    else
        [s,mess,messid]=mkdir(directorio);
        if s==1
            clc
            disp('   ')
            disp('Se creo el directorio: ')
            disp(directorio)
        else
            disp('   ');
            disp(mess)
            dips('   ')
            dips(messid)
        end
    end
end


