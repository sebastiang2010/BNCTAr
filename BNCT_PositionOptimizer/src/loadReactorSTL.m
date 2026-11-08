function [wallStruct, bodyStruct] = loadReactorSTL(wallPath, bodyPath)
%LOADREACTORSTL Carga la geometria del reactor (puerto del beamline).
%   [wallStruct, bodyStruct] = LOADREACTORSTL(WALLPATH, BODYPATH)
%
%   WALLPATH : malla de la pared/estructura cercana al paciente
%              (formatos soportados: .stl nativo via stlread; .vtk ASCII).
%   BODYPATH : (OPCIONAL) malla adicional del cuerpo del reactor. Si se
%              omite o esta vacia, bodyStruct se devuelve vacio y el SDF /
%              colisiones usan solo wallStruct.
%
%   SALIDA NORMALIZADA (ambas):
%     .points : Nx3 double (mm)
%     .faces  : Mx3 double, indices de triangulo 1-based (referencia a .points)
%     .fuente : ruta/archivo de origen (informativo)
%
%   NOTAS:
%     - MATLAB stlread devuelve un objeto triangulation (NO un struct con
%       vertices/faces): esta funcion normaliza ambos casos.
%     - .vtk: se implementa un parser ASCII minimo (POINTS + POLYGONS/CELLS
%       con celdas de 3 o 4 vertices; cuadrilateros se dividen en 2
%       triangulos). VTK binario o .vtp (XML) -> error claro.
%     - Unidades: se asume que la malla esta en mm (marco del haz).

if nargin < 2
    bodyPath = '';
end

wallStruct = cargarMalla(wallPath, 'wall');

if isempty(bodyPath)
    bodyStruct = struct('points', zeros(0,3), 'faces', zeros(0,3), ...
        'fuente', '(no provisto)');
else
    bodyStruct = cargarMalla(bodyPath, 'body');
end
end

% ------------------------------------------------------------------------
function outStruct = cargarMalla(path, etiqueta)
% Carga y normaliza una malla segun su extension.
if isempty(path) || exist(path, 'file') ~= 2
    error('BNCT:loadReactorSTL:archivoNoEncontrado', ...
        'Malla "%s" no encontrada: %s', etiqueta, path);
end
[~,~,ext] = fileparts(path);
switch lower(ext)
    case '.stl'
        outStruct = normalizar(stlread(path));
    case '.vtk'
        outStruct = parseVTKascii(path);
    otherwise
        error('BNCT:loadReactorSTL:formatoNoSoportado', ...
            ['Extension "%s" no soportada para "%s". Use .stl (nativo) o ' ...
             '.vtk ASCII (POINTS + POLYGONS/CELLS).'], ext, etiqueta);
end
outStruct.fuente = path;
end

% ------------------------------------------------------------------------
function outStruct = normalizar(mesh)
% Normaliza triangulation o struct (vertices/faces) a .points/.faces 1-based.
if isa(mesh, 'triangulation')
    outStruct.points = double(mesh.Points);
    outStruct.faces  = double(mesh.ConnectivityList);
elseif isstruct(mesh) && (isfield(mesh, 'vertices') || isfield(mesh, 'points'))
    if isfield(mesh, 'points')
        outStruct.points = double(mesh.points);
    else
        outStruct.points = double(mesh.vertices);
    end
    if isfield(mesh, 'faces')
        outStruct.faces = double(mesh.faces);
    elseif isfield(mesh, 'ConnectivityList')
        outStruct.faces = double(mesh.ConnectivityList);
    else
        error('BNCT:loadReactorSTL:structSinCaras', ...
            'El struct de malla no tiene campo "faces" ni "ConnectivityList".');
    end
else
    error('BNCT:loadReactorSTL:tipoNoReconocido', ...
        ['stlread devolvio un tipo no reconocido (%s). Se esperaba ' ...
         'triangulation o struct.'], class(mesh));
end
validarMalla(outStruct);
end

% ------------------------------------------------------------------------
function outStruct = parseVTKascii(path)
% Parser ASCII minimo para VTK legacy POLYDATA (POINTS + POLYGONS/CELLS).
% Lee cada seccion acotado por el conteo declarado en su encabezado y
% IGNORA cualquier otra seccion (LINES, VERTICES, POINT_DATA, FIELD...),
% de modo que los datos nunca se contaminan entre secciones.
fid = fopen(path, 'rt');
if fid == -1
    error('BNCT:loadReactorSTL:noSePuedeAbrir', 'No se pudo abrir: %s', path);
end
texto = fread(fid, '*char')';
fclose(fid);

if contains(texto, 'BINARY')
    error('BNCT:loadReactorSTL:vtkBinario', ...
        ['"%s" es VTK BINARIO. Este parser solo lee ASCII. Convierta el ' ...
         'archivo a STL o a VTK ASCII (p.ej. con meshlab/paraview).'], path);
end

lineas = regexp(texto, '\r?\n', 'split');
nLin = numel(lineas);
patronFloat = '[-+]?\d+\.?\d*(?:[eEdD][-+]?\d+)?';

points = zeros(0, 3);
faces  = zeros(0, 3);
i = 0;
while i < nLin
    i = i + 1;
    L = strtrim(lineas{i});
    tok = regexp(L, '^([A-Za-z_]\w*)\s+(\d+)\s*(\w*)', 'tokens', 'once');
    if isempty(tok)
        continue;                          % linea de datos o comentario
    end
    switch upper(tok{1})
        case 'POINTS'
            nPuntos = str2double(tok{2});
            [vals, i] = leerNumeros(lineas, i + 1, 3 * nPuntos, patronFloat);
            if numel(vals) < 3 * nPuntos
                error('BNCT:loadReactorSTL:vtkPuntosIncompletos', ...
                    'POINTS declara %d puntos pero el archivo tiene menos numeros.', nPuntos);
            end
            points = reshape(vals(1:3*nPuntos), 3, nPuntos)';   % Nx3

        case {'POLYGONS', 'CELLS'}
            numCeldas = str2double(tok{2});
            nDatos    = str2double(tok{3});
            [vals, i] = leerNumeros(lineas, i + 1, nDatos, '\d+');
            fcsTmp = zeros(2 * numCeldas, 3);  % max 2 triangulos por celda
            nTr = 0;
            k = 1;
            while k <= numel(vals) && nTr < 2 * numCeldas
                nv = vals(k);
                if nv < 1 || k + nv > numel(vals)
                    break;                     % datos truncos: descartar cola
                end
                idx = vals(k+1 : k+nv) + 1;    % VTK 0-based -> 1-based
                if idx(end) > size(points, 1) || any(idx < 1)
                    warning('BNCT:loadReactorSTL:indiceInvalido', ...
                        'Celda con indices fuera de rango en "%s" (descartada).', path);
                elseif nv == 3
                    nTr = nTr + 1;
                    fcsTmp(nTr, :) = idx(:).';
                elseif nv == 4
                    nTr = nTr + 1;
                    fcsTmp(nTr, :) = idx(1:3);
                    nTr = nTr + 1;
                    fcsTmp(nTr, :) = idx([1 3 4]);    % quad -> 2 triangulos
                else
                    warning('BNCT:loadReactorSTL:celdaIgnorada', ...
                        'Celda de %d vertices ignorada (solo 3 o 4).', nv);
                end
                k = k + nv + 1;
            end
            faces = fcsTmp(1:nTr, :);

        otherwise
            % Seccion no soportada (LINES, VERTICES, POINT_DATA, ...):
            % se ignora por completo; su contenido se procesa cuando su
            % propio encabezado aparezca (nunca, para las no soportadas).
            continue;
    end
end

if isempty(points)
    error('BNCT:loadReactorSTL:vtkSinPoints', ...
        'No se encontro la seccion POINTS en "%s".', path);
end
if isempty(faces)
    error('BNCT:loadReactorSTL:vtkSinTriangulos', ...
        ['No se extrajo ningun triangulo de "%s" (POLYGONS/CELLS ' ...
         'ausente o sin celdas de 3-4 vertices).'], path);
end
outStruct = struct('points', points, 'faces', faces);
validarMalla(outStruct);
end

% ------------------------------------------------------------------------
function [vals, i] = leerNumeros(lineas, i0, cuantos, patron)
% Acumula exactamente CUANTOS numeros leidos desde la linea i0. Se detiene
% al completar el conteo o al encontrar la proxima linea-encabezado.
vals = zeros(0, 1);
i = i0 - 1;
nLin = numel(lineas);
while i < nLin && numel(vals) < cuantos
    i = i + 1;
    L = strtrim(lineas{i});
    if isempty(L) || L(1) == '#'
        continue;
    end
    % un encabezado de otra seccion corta la lectura
    if ~isempty(regexp(L, '^[A-Za-z_]\w*\s+\d+', 'once')) && isempty(regexp(L, '^\d', 'once'))
        i = i - 1;                         % devolver la linea para el bucle principal
        break;
    end
    trozos = regexp(L, patron, 'match');
    trozos = regexprep(trozos, '[dD]', 'E');   % exponente Fortran -> estandar
    vals = [vals; str2double(trozos(:))]; %#ok<AGROW>
end
end

% ------------------------------------------------------------------------
function validarMalla(S)
if size(S.points, 2) ~= 3
    error('BNCT:loadReactorSTL:puntosMalFormados', ...
        'points debe ser Nx3; recibido %dx%d.', size(S.points,1), size(S.points,2));
end
if size(S.faces, 2) ~= 3
    error('BNCT:loadReactorSTL:carasMalFormadas', ...
        'faces debe ser Mx3; recibido %dx%d.', size(S.faces,1), size(S.faces,2));
end
if ~isempty(S.faces) && (max(S.faces(:)) > size(S.points,1) || min(S.faces(:)) < 1)
    error('BNCT:loadReactorSTL:indicesFueraDeRango', ...
        'faces referencia indices fuera de [1..N] de points.');
end
end
