function [D1c, meta] = leerPHITSGamma(rutaArchivo)
%LEERPHITSGAMMA Lee un .out PHITS de fotones (1 mset, T-Deposit o T-Track).
%   [D1c, meta] = LEERPHITSGAMMA(RUTA_ARCHIVO)
%
%   Adaptado NO-INTERACTIVO de f_lectura_phits_gamma.m del directorio padre
%   (codigo del usuario, AUTORIDAD DE SEMANTICA). Se copia su logica casi
%   literal; las unicas diferencias estan listadas abajo y marcadas en el
%   codigo con "ADAPTACION"/"CORRECCION".
%
%   Secciones copiadas de f_lectura_phits_gamma.m (159 lineas):
%     - lineas 33-45 : tokens (' newpage:', ' part ', ' nx =', ...) y part.
%     - lineas 48-63 : nx, ny, nz via regexp.
%     - lineas 69-79 : deteccion de paginas por lineas 'hc' (a1(1:2:end),
%                      porque cada pagina trae DOS lineas con 'hc': el
%                      descriptor de datos y el del pie) y colum (lineas
%                      por pagina).
%     - lineas 114-143: lectura de cada pagina y armado C(:,:,k).
%     - lineas 152-154: D1c.Dg = C si part es photon.
%
%   ADICION (marcada): el original NO parsea la geometria (xmin..zdel) y
%   no devuelve xc/yc/zc. Para que loadPHITSDose pueda armar la grilla en
%   mm se copia el bloque de geometria de f_lectura_phits_neutron.m
%   (lineas 38-79), que es el mismo formato de header PHITS.
%
%   DIFERENCIAS RESPECTO AL ORIGINAL (obligadas por los datos reales,
%   grilla nx=100, ny=200, nz=200, 200 paginas):
%     A1) RUTA_ARCHIVO entra como argumento en lugar de uigetfile().
%     A2) Devuelve (D1c, meta); elimina clc/cd/figuras interactivas.
%     F1) ADAPTACION DE RENDIMIENTO: conversion texto->numero por bloque
%         con sscanf('%f'), con fallback al camino literal (cellfun de
%         strsplit+str2double con padding NaN) si el bloque trae tokens
%         no numericos. Equivalente numericamente (verificado).
%     F2) CORRECCION de conteo de lineas: el original lee
%         lines(a(k)+1 : a(k)+colum-1), que con grillas donde
%         mod(nx*ny,10)==0 lee UNA LINEA DE MENOS y el reshape posterior
%         tira error (verificado en la grilla real 100x200: 19990 vs
%         20000 valores). La cantidad correcta es ceil(nx*ny/10).
%     F3) CORRECCION de la transposicion: el original hace
%         C(:,:,k)=B1' con B1=B' en la rama divisible (neto: correcto)
%         pero B1=B en la rama no divisible (neto: transpuesto, y con
%         nx~=ny redimensiona C implicitamente a ny x nx). Se unifica:
%         la pagina se arma x-rapida (B1 = reshape/loop, igual que
%         leerPHITSNeutron) y se guarda C(:,:,k)=B1.
%     F4) CORRECCION de orientacion en y (fliplr): el archivo declara
%         "# ( ( data(x,y), x = 1, nx ), y = ny, 1, -1 )": y DESCiende
%         en el stream mientras yc asciende. Igual correccion que
%         bnct.py (flipud) y que leerPHITSNeutron.
%
%   SALIDAS
%     D1c : struct con
%       .Dg : [nx ny nz] Gy/source (fotones)
%       .xc/.yc/.zc : centros de voxel en cm (adicion documentada arriba)
%     meta : struct con .part .nx .ny .nz .xmin .xmax .ymin .ymax .zmin
%            .zmax .xdel .ydel .zdel .nPaginas .archivo .tiempoParseo

tInicio = tic;

% ======================================================================
% (lineas 33-63 del original) tokens, part, nx/ny/nz
% ======================================================================
if nargin < 1 || isempty(rutaArchivo)
    error('BNCT:leerPHITSGamma:faltaRuta', 'Se requiere la ruta del archivo .out (modo no interactivo).');
end
if exist(rutaArchivo, 'file') ~= 2
    error('BNCT:leerPHITSGamma:noExiste', 'No existe el archivo: %s', rutaArchivo);
end

lines=readlines(rutaArchivo);

txt=' newpage:';
txt1=' part ';
txt2=' nx =';
txt3=' ny =';
txt4=' nz =';

% busco particula
b=find(contains(lines,txt1));
p = strtrim(lines(b(1)));
p = strsplit(p, '=');
% Divide por el signo '='
part=strtrim(p{2});
clear p b %#ok<NASGU>

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

clear c p %#ok<NASGU>

% ======================================================================
% ADICION (de f_lectura_phits_neutron.m, lineas 38-79): geometria
% ======================================================================
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
% (lineas 69-79 del original) paginas por lineas 'hc'
% ======================================================================
a1=find(contains(lines,'hc'));
a=a1(1:2:end); % porque esta entre medio hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;

% CORRECCION F2: cantidad de lineas de datos por pagina = ceil(nx*ny/10).
% (El original usaba colum = round(nx*ny/10) (+1 si hay resto) y leia
%  colum-1 lineas: correcto para grillas no divisibles por 10, pero
%  lea una linea de menos en grillas divisibles -> error de reshape.)
colum = ceil(nx*ny/10); %#ok<NASGU>  (se conserva el nombre del original)

if numel(a) < nz
    error('BNCT:leerPHITSGamma:paginasInsuficientes', ...
        'Se esperaban %d paginas y el archivo tiene %d marcadores hc utiles.', nz, numel(a));
end

% ======================================================================
% (lineas 114-143 del original) lectura de cada pagina
% ======================================================================
C=ones(nx,ny,nz);
A=[];
for k=1:nz
    % CORRECCION F2: leer 'colum' lineas (el original leia colum-1)
    A{k}=lines(a(k)+1:a(k)+colum);

    % F1) conversion por bloque (sscanf) con fallback literal
    valores = convertirBloque(A{k}, colum, nx*ny);

    % F3/F4) pagina x-rapida + volteo en y (ver cabecera)
    B1 = armarPagina(valores, nx, ny);
    B1 = fliplr(B1);

    C(:,:,k)=B1;
end

% ======================================================================
% (lineas 152-154 del original) empaquetado
% ======================================================================
if ~strcmpi(part, 'photon')
    warning('BNCT:leerPHITSGamma:partInesperada', ...
        'El archivo declara part=%s (se esperaba photon). Revisar el tally.', part);
end
D1c.Dg=C;

D1c.xc=xc;
D1c.yc=yc;
D1c.zc=zc;

meta = struct('part', part, 'nx', nx, 'ny', ny, 'nz', nz, ...
    'xmin', xmin, 'xmax', xmax, 'ymin', ymin, 'ymax', ymax, ...
    'zmin', zmin, 'zmax', zmax, 'xdel', xstep, 'ydel', ystep, 'zdel', zstep, ...
    'nPaginas', numel(a), 'archivo', char(rutaArchivo), ...
    'tiempoParseo', toc(tInicio));
end

% ========================================================================
function valores = convertirBloque(A, mLineas, nEsperado)
% F1) Conversion texto->numero del bloque de la pagina (misma estrategia
%     que leerPHITSNeutron): sscanf rapido + fallback literal del usuario
%     (cellfun strsplit+str2double, padding NaN de la ultima linea corta,
%     identico a f_lectura_phits_gamma.m lineas 116-121).
bloque = strjoin(string(A), newline);
valores = sscanf(bloque, '%f');
valores = valores(:).';
if numel(valores) ~= nEsperado
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
function B1 = armarPagina(valores, nx, ny)
% Pagina 2D [nx x ny] en orden x-rapido (semantica de las dos ramas del
% original unificada, CORRECCION F3):
%   - divisible: reshape(valores, nx, ny) == reshape(B13',nx,ny) del
%     original, que ya es x-rapido.
%   - no divisible: bucle B(:,j)=valores[(j-1)*nx+(1:nx)], igual que las
%     lineas 131-136 del original.
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
