function MeshData =lecturaDosis(filename,dir)
% *************************************************************************
% ABRO EL GUI DE OPENFILE Y CARGO EL ARCHIVO. SI EL USUARIO PRESIONA CANCEL
% DEVUELVO UNA ESTRUCTURA NULA
% *************************************************************************
switch(nargin)
    case 0
        [filename_1,directorio]=uigetfile('*.*','Enter MCNP neutron output file');
    case 2
        filename_1=filename;
        directorio=dir;
end
if isequal(filename_1,0) || isequal(directorio,0)
    disp('User pressed cancel')
    MeshData=0;
else
    path(path,directorio);
    file1= fopen(filename_1, 'rt');
    % [FILENAME,PERMISSION,MACHINEFORMAT,ENCODING]=fopen(file1)
    % *************************************************************************
    % INICIALIZO LAS VARIABLES A USAR
    % *************************************************************************
    y=0;
    MeshData.X0=0;
    MeshData.Y0=0;
    MeshData.Z0=0;
    MeshData.Xf=0;
    MeshData.Yf=0;
    MeshData.Zf=0;
    MeshData.StepX=0;
    MeshData.StepY=0;
    MeshData.StepZ=0;
    MeshData.FlagBoro=0;
    MeshData.FlagTermico=0;
    MeshData.FlagRapido=0;
    MeshData.FlagGamma=0;
    MeshData.FlagFlujo=0;
    fotones=0;
    % *************************************************************************
    % FINALIZO INICIALIZACION DE VARIABLES - COMIENZO A RECORRER EL ARCHIVO
    % *************************************************************************
    while feof(file1) == 0 %indica 1 cuando es el final del archivo
        tline = fgetl(file1);%lee linea por linea
        % *************************************************************************
        % SI LA MESH ES DE PHOTONS COLOCO LA VARIABLE fotones EN 1
        % *************************************************************************
        % if (findstr(tline,' This is a photon mesh tally.'))
        %     fotones=1;
        % elseif (findstr(tline,' photon   mesh tally.'))
        %     fotones=1;
        % end
        % Reemplaza múltiples espacios por uno solo
        tline_clean = regexprep(tline, '\s+', ' ');

        if contains(tline_clean, 'This is a photon mesh tally.') || ...
                contains(tline_clean, 'photon mesh tally.')
            fotones = 1;
        end

        % *************************************************************************
        % BUSCO EL ENCABEZADO QUE INDICA INICIO DE LOS DATOS EN EL ARCHIVO
        % *************************************************************************
        %matches = strfind(tline, '        X         Y         Z     Result     Rel Error');
        matches=regexp(tline, 'X\s+Y\s+Z\s+Result\s+Rel\s+Error', 'once');
        if ~isempty(matches)
            %if matches > 0 % si encuntra devuelve la posicion>0
            if(strfind(tline,'Energy'))
                Energy=1;

            else
                Energy=0;
            end

            a=fscanf(file1,'%g',[5+Energy,inf]); %lee en el archivo find1,
            a=a';
            d1=0;
            d2=0;
            MeshData.X0=a(1,1+Energy);
            MeshData.Y0=a(1,2+Energy);
            MeshData.Z0=a(1,3+Energy);

            L=length(a(:,1));
            if (Energy==1)
                Bins=a(:,1)-a(1,1); % Solo es cero cuando se repite la secuencia
                Nbins=numel(unique(Bins));
                fin=L/Nbins;
            else
                fin=L;
                Nbins=1;
            end
            truco1=a(1:fin,1+Energy)-a(1,1+Energy); % Solo es cero cuando se repite la secuencia
            truco2=a(1:fin,2+Energy)-a(1,2+Energy); %Solo es cero cuando se repite la secuencia en Y

            d1=numel(unique(a(1:fin,1+Energy)));%sum(truco1==0);
            Xu=unique(a(1:fin,1+Energy));
            d2=numel(unique(a(1:fin,2+Energy))); %(sum(truco2==0);
            Yu=unique(a(1:fin,2+Energy));
            %                i=1;
            %                d2=0
            %                while i~=-1
            %                    if truco2(i,1)==0
            %                      d2=d2+1; %Cuenta las veces que se repite la secuencia
            %                      i=i+1;
            %                    else
            %                      MeshData.StepY=truco2(i,1);
            %                      i=-1;
            %                    end;
            %                end;
            MeshData.StepX=Xu(2)-Xu(1);
            MeshData.StepY=Yu(2)-Yu(1);
            MeshData.StepZ=(a(2,3+Energy)-a(1,3+Energy));
            MeshData.tamX=d1;
            MeshData.tamY=d2;
            MeshData.tamZ=fin/MeshData.tamX/MeshData.tamY;
            % *************************************************************************
            % ASIGNO CADA COMPONENTE DE DOSIS A SU RESPECTIVA VARIABLE
            % *************************************************************************
            % BORO CORRESPONDE A DOSIS BORO

            for idx=1:Nbins
                y=y+1;
                bgn=(idx-1)*(L/Nbins)+1;
                fsh=idx*(L/Nbins);

                if (fotones==0)
                    Data=permute(reshape(a(bgn:fsh,4+Energy),MeshData.tamZ,MeshData.tamY,MeshData.tamX),[3 2 1]);
                    Error=permute(reshape(a(bgn:fsh,5+Energy),MeshData.tamZ,MeshData.tamY,MeshData.tamX),[3 2 1]);
                    MeshData.type(y)=cellstr('n ');

                else
                    Data=permute(reshape(a(bgn:fsh,4+Energy),MeshData.tamZ,MeshData.tamY,MeshData.tamX),[3 2 1]);
                    Error=permute(reshape(a(bgn:fsh,5+Energy),MeshData.tamZ,MeshData.tamY,MeshData.tamX),[3 2 1]);
                    MeshData.type(y)=cellstr('p ');
                end
                MeshData.Matrix(:,:,:,y)=Data;
                MeshData.Error(:,:,:,y)=Error;

                if(Energy==1)
                    MeshData.En(y)=a(bgn,1);
                else
                    MeshData.En(y)=inf;
                end
            end




        end
    end
    MeshData.Comps=y;
    % *************************************************************************
    %EXTRAIGO DATOS DE DIMENSIONES DE LA MESH
    % *************************************************************************


    MeshData.Xf=(MeshData.tamX*MeshData.StepX+MeshData.X0)-MeshData.StepX;
    MeshData.Yf=(MeshData.tamY*MeshData.StepY+MeshData.Y0)-MeshData.StepY;
    MeshData.Zf=(MeshData.tamZ*MeshData.StepZ+MeshData.Z0)-MeshData.StepZ;
    MeshData.ejeX=(MeshData.X0:MeshData.StepX:MeshData.Xf);
    MeshData.ejeY=(MeshData.Y0:MeshData.StepY:MeshData.Yf);
    MeshData.ejeZ=(MeshData.Z0:MeshData.StepZ:MeshData.Zf);
    MeshData.Ro=[MeshData.X0 MeshData.Y0 MeshData.Z0];
    MeshData.Rf=[MeshData.Xf MeshData.Yf MeshData.Zf];
    MeshData.Step=[MeshData.StepX MeshData.StepY MeshData.StepZ];
    fclose(file1);
end





