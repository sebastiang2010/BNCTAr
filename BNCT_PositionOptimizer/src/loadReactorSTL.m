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

% --- puntos ---
tok = regexp(texto, 'POINTS\s+(\d+)\s+(\w+)', 'tokens', 'once');
if isempty(tok)
    error('BNCT:loadReactorSTL:vtkSinPoints', ...
        'No se encontro la seccion POINTS en "%s".', path);
end
nPuntos = str2double(tok{1});
% los 3*nPuntos numeros siguen al encabezado POINTS: los localizamos por
% posicion (buscamos el indice del token "POINTS" en la lista cruda).
posPOINT = strfind(texto, 'POINTS');
cola = texto(posPOINT(1):end);
todos = regexp(cola, '[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?', 'match');
% los primeros 2 tokens son "n" y "tipo"; luego 3*nPuntos numeros
necesarios = 2 + 3*nPuntos;
if numel(todos) < necesarios
    error('BNCT:loadReactorSTL:vtkPuntosIncompletos', ...
        'POINTS declara %d puntos pero el archivo tiene menos numeros.', nPuntos);
end
vals = str2double(todos(3:2+3*nPuntos));
points = reshape(vals, 3, nPuntos)';      % Nx3

% --- poligonos ---
tok = regexp(texto, 'POLYGONS\s+(\d+)\s+(\d+)', 'tokens', 'once');
if isempty(tok)
    tok = regexp(texto, 'CELLS\s+(\d+)\s+(\d+)', 'tokens', 'once');
end
if isempty(tok)
    error('BNCT:loadReactorSTL:vtkSinPoligonos', ...
        ['No se encontro la seccion POLYGONS/CELLS en "%s". Si la malla ' ...
         'contiene solo LINES/STRIPS no es usable como superficie.'], path);
end
posPOLY = strfind(texto, 'POLYGONS');
if isempty(posPOLY); posPOLY = strfind(texto, 'CELLS'); end
cola = texto(posPOLY(1):end);
todos = regexp(cola, '\d+', 'match');
necesarios = 2 + str2double(tok{2});
if numel(todos) < necesarios
    error('BNCT:loadReactorSTL:vtkPoligonosIncompletos', ...
        'POLYGONS/CELLS declara mas datos de los presentes en el archivo.');
end
vals = str2double(todos(3:2+str2double(tok{2})));

% recorrer celdas: primera cifra = nVertices de la celda
numCeldas = str2double(tok{1});
fcsTmp = zeros(2 * numCeldas, 3);   % preallocado: max 2 triangulos por celda
nTr = 0;
i = 1;
while i <= numel(vals)
    nv = vals(i);
    if i + nv > numel(vals) + 1
        break; % datos truncos: se descarta la celda incompleta
    end
    idx = vals(i+1 : i+nv) + 1;           % VTK es 0-based -> 1-based
    if idx(end) > nPuntos || any(idx < 1)
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
    i = i + nv + 1;
end
faces = fcsTmp(1:nTr, :);

if isempty(faces)
    error('BNCT:loadReactorSTL:vtkSinTriangulos', ...
        'No se extrajo ningun triangulo de "%s".', path);
end
outStruct = struct('points', points, 'faces', faces);
validarMalla(outStruct);
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
