function [D1c, meta] = leerPHITSNeutron(rutaArchivo)
%LEERPHITSNEUTRON Lee un .out PHITS [T-Track] de neutrones (3 msets).
%   [D1c, meta] = LEERPHITSNEUTRON(RUTA_ARCHIVO)
%
%   Adaptado NO-INTERACTIVO de f_lectura_phits_neutron.m del directorio
%   padre (codigo del usuario, AUTORIDAD DE SEMANTICA). Se copia su logica
%   casi literal; las unicas diferencias se listan abajo y estan marcadas
%   con "ADAPTACION"/"CORRECCION" en el codigo.
%
%   Secciones copiadas de f_lectura_phits_neutron.m (348 lineas):
%     - lineas 27-79  : readlines + tokens de geometria (xmin..zdel) y
%                       centros de voxel xc/yc/zc (queda en cm, igual que
%                       el original; la conversion a mm la hace el
%                       llamador, p.ej. loadPHITSDose).
%     - lineas 84-116 : deteccion de part y lectura de nx, ny, nz.
%     - lineas 118-134: indices de '#newpage:'/' newpage:', m=round(nx*ny/10),
%                       salto=25 lineas de cabecera por pagina.
%     - lineas 136-276: bucles de los 3 msets (1=Boron -> .boro,
%                       2=FastNeutron -> .fast, 3=ThermalNeutron -> .thn)
%                       en ORDEN DE APARICION en el archivo (verificado:
%                       200 paginas por mset, bloques contiguos).
%     - lineas 326-345: empaquetado D1c (.boro/.fast/.thn + .xc/.yc/.zc).
%
%   DIFERENCIAS RESPECTO AL ORIGINAL (todas obligadas por los datos reales
%   de BNCT_PositionOptimizer/data/, grilla nx=100, ny=200, nz=200):
%     A1) RUTA_ARCHIVO entra como argumento en lugar de uigetfile()
%         (necesario para automatizacion).
%     A2) Todo se devuelve como salidas (D1c, meta); se eliminan clc/cd/
%         disp interactivos. Las figuras/comentarios del original no
%         aplican en modo automatico.
%     F1) ADAPTACION DE RENDIMIENTO: la conversion texto->numero se hace
%         por bloque con sscanf('%f') y, si el bloque trae tokens no
%         numericos (p.ej. '*****'), cae al camino LITERAL del usuario
%         (cellfun de strsplit+str2double con padding NaN de la ultima
%         linea). Con cellfun literal, las 600 paginas tardan ~300 s;
%         con sscanf, ~15 s. Verificado igualdad numerica en paginas de
%         prueba del archivo real.
%     F2) CORRECCION de la rama divisible (mod(nx*ny,10)==0): el original
%         hace B=reshape(B13',nx,ny); B1=B'; y luego C(:,:,n,g)=B1. Con
%         nx~=ny eso CRASHEA (B1 es ny x nx y el slice es nx x ny; en la
%         grilla real 100x200 tira "Subscripted assignment dimension
%         mismatch"), y ademas B1=B' contradice la propia rama no
%         divisible (que deja B1=B). reshape(B13',nx,ny) YA produce el
%         orden x-rapido correcto (verificado con ejemplo minimo y contra
%         el parser bnct.py del usuario). Se usa B1=reshape(B13',nx,ny).
%     F3) CORRECCION de orientacion en y (fliplr): el propio archivo
%         declara "# ( ( data(x,y), x = 1, nx ), y = ny, 1, -1 )" e
%         imprime por pagina "hc:  y = +29.85 to -29.85 by 0.3 ; x = ...",
%         es decir, y DESCiende a lo largo del stream mientras xc/yc/zc
%         ascienden. El original deja B1 sin voltear, lo que espeja la
%         dosis en y respecto de sus propias etiquetas yc. bnct.py
%         (parsePHITS, traduccion vigente en Slicer) corrige lo mismo con
%         flipud(B.T), que es identico a fliplr(B) transpuesto. Se aplica
%         fliplr para que D1c.boro(i,j,k) sea la dosis en (xc(i),yc(j),zc(k)).
%
%   SALIDAS
%     D1c : struct con
%       .boro : [nx ny nz] Gy/source, mset 1  (Boron)
%       .fast : [nx ny nz] Gy/source, mset 2  (FastNeutron)
%       .thn  : [nx ny nz] Gy/source, mset 3  (ThermalNeutron)
%       .xc/.yc/.zc : centros de voxel en cm (convencion del original)
%     meta : struct con .part .nx .ny .nz .xmin .xmax .ymin .ymax .zmin
%            .zmax .xdel .ydel .zdel .nPaginas .archivo .tiempoParseo
%
%   Ver README.md (seccion 6) para el uso desde loadPHITSDose.

tInicio = tic;

% ======================================================================
% (lineas 27-79 del original) lectura y geometria
% ======================================================================
if nargin < 1 || isempty(rutaArchivo)
    error('BNCT:leerPHITSNeutron:faltaRuta', 'Se requiere la ruta del archivo .out (modo no interactivo).');
end
if exist(rutaArchivo, 'file') ~= 2
    error('BNCT:leerPHITSNeutron:noExiste', 'No existe el archivo: %s', rutaArchivo);
end

lines = readlines(rutaArchivo);   % igual que el original

txt=' newpage:';
txt1=' part ';
txt2=' nx =';
txt3=' ny =';
txt4=' nz =';
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
    datos_mesh{1,i} = str2double(num{1}); %#ok<AGROW>
end
xmin = datos_mesh{1,1};
xmax = datos_mesh{1,2};
ymin = datos_mesh{1,4};
ymax = datos_mesh{1,5};
zmin = datos_mesh{1,7};
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

% ======================================================================
% (lineas 84-116 del original) particula y nx/ny/nz
% ======================================================================
b=find(contains(lines,txt1));
p = strtrim(lines(b(1)));
p = strsplit(p, '=');
part=strtrim(p{2});
clear p b %#ok<NASGU>

% Aviso de seguridad en modo automatico (el original solo avisaba en
% pantalla; si el archivo no es de neutrones, los msets no significan nada)
if ~strcmpi(part, 'neutron')
    warning('BNCT:leerPHITSNeutron:partInesperada', ...
        'El archivo declara part=%s (se esperaba neutron). Revisar el tally.', part);
end

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

clear c p %#ok<NASGU>

% ======================================================================
% (lineas 118-134 del original) indices de paginas y constantes
% ======================================================================
a=find(contains(lines,txt));
a1=find(lines=='#newpage:'); % la primera es diferente

a=[a1,a'];

nPaginas = numel(a);
if nPaginas < 3*nz
    error('BNCT:leerPHITSNeutron:paginasInsuficientes', ...
        ['Se esperaban %d paginas (3 msets x %d slices en z) y el archivo ' ...
         'tiene %d marcadores newpage.'], 3*nz, nz, nPaginas);
end

m=round(nx*ny/10);
salto=25; %desde newpage a donde empieza la matriz (verificado: 'hc:' esta en
          % newpage+24 y los datos arrancan en newpage+25 en las 600 paginas)
m1=m-1; % tamano de la lectura

% ======================================================================
% (lineas 136-276 del original) bucles de los 3 msets
% ======================================================================
C=ones(nx,ny,nz,3);

for g=1:3

    nmset=g;
    switch g
        case 1
            rangoK = 1:nz;              % paginas 1..nz      -> mset 1
        case 2
            rangoK = nz+1:2*nz;         % paginas nz+1..2nz -> mset 2
        case 3
            rangoK = 2*nz+1:3*nz;       % paginas 2nz+1..3nz-> mset 3
    end

    n=1;
    for k=rangoK

        A=lines(a(k)+salto:a(k)+salto+m1);

        % F1) conversion por bloque (sscanf) con fallback literal
        valores = convertirBloque(A, m, nx*ny);

        % F2/F3) armado de la pagina (ver comentarios de cabecera)
        B = armarPagina(valores, nx, ny, m);
        B1 = fliplr(B);   % F3: el archivo escribe y descendente

        C(:,:,n,nmset)=B1;
        n=n+1;
    end
end

% ======================================================================
% (lineas 326-345 del original) empaquetado de salida
% ======================================================================
for i=1:3
    A=C(:,:,:,i);
    A1=A;
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

meta = struct('part', part, 'nx', nx, 'ny', ny, 'nz', nz, ...
    'xmin', xmin, 'xmax', xmax, 'ymin', ymin, 'ymax', ymax, ...
    'zmin', zmin, 'zmax', zmax, 'xdel', xstep, 'ydel', ystep, 'zdel', zstep, ...
    'nPaginas', nPaginas, 'archivo', char(rutaArchivo), ...
    'tiempoParseo', toc(tInicio));
end

% ========================================================================
function valores = convertirBloque(A, mLineas, nEsperado)
% F1) Conversion texto->numero de las mLineas de una pagina.
%     Camino rapido: un sscanf('%f') sobre el bloque concatenado (los
%     datos PHITS son todos numericos en notacion E). Si la cantidad no
%     da exacta (token no numerico, p.ej. '*****'), cae al camino LITERAL
%     del usuario: cellfun de strsplit+str2double por linea con padding
%     NaN de la ultima linea corta (identico a f_lectura_phits_neutron.m
%     lineas 157-165).
bloque = strjoin(string(A), newline);
valores = sscanf(bloque, '%f');
valores = valores(:).';   % fila, igual que el orden del stream
if numel(valores) ~= nEsperado
    % camino literal del usuario (lento, pero tolera tokens raros)
    B10 = cell2mat(cellfun(@(x) str2double(strsplit(strtrim(x))), A(1:end-1), 'UniformOutput', false));
    B11 = cell2mat(cellfun(@(x) str2double(strsplit(strtrim(x))), A(end,:), 'UniformOutput', false));
    B12 = NaN(1,10);
    B12(1:length(B11))=B11;
    B13=[B10;B12];
    valores = reshape(B13.', 1, []);
    valores = valores(1:min(nEsperado, numel(valores)));
    if numel(valores) < nEsperado
        valores(end+1:nEsperado) = NaN; %#ok<AGROW>
    end
end
end

% ========================================================================
function B1 = armarPagina(valores, nx, ny, mLineas)
% Armado de la pagina 2D [nx x ny] a partir del stream de mLineas x 10
% valores, siguiendo la semantica del original:
%   - rama no divisible (lineas 171-181): B(:,j) toma los nx valores
%     siguientes  ->  B(i,j) = valor[(j-1)*nx + i]  (x mas rapido).
%   - rama divisible (lineas 167-169): reshape(B13',nx,ny), que produce
%     exactamente la misma matriz (verificado). Se corrige B1=B' -> B1=B
%     (CORRECCION F2) para que ambas ramas coincidan y no crashee con
%     nx ~= ny en C(:,:,n,g)=B1.
% El padding de la ultima linea corta ya lo hizo convertirBloque (NaN).
if mod(nx*ny,10)==0
    B1 = reshape(valores, nx, ny);
else
    B = zeros(nx, ny);
    for j=1:ny
        B(:,j) = valores((j-1)*nx + (1:nx));
    end
    B1 = B;
end
end
