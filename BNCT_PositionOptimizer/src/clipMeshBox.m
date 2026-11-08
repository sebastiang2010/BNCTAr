function [S2, kept] = clipMeshBox(S, box)
% CLIPMESHBOX Conserva solo los triangulos completamente dentro de una caja.
%
%   [S2, kept] = clipMeshBox(S, box)
%
% ENTRADAS
%   S    - struct de malla con campos .points (Nx3, mm) y .faces (Mx3,
%          indices 1-based). Se acepta .triangles como alias de .faces.
%   box  - [xmin xmax ymin ymax zmin zmax] en mm.
%
% SALIDAS
%   S2   - malla recortada (.points/.faces); los puntos que quedaron sin
%          referencia se descartan y los indices se renumera.
%   kept - indices (M2x1) de los triangulos originales conservados.
%
% CRITERIO
%   Un triangulo se conserva solo si sus TRES vertices estan dentro de la
%   caja. No se recortan triangulos parciales: en mallas de colision un
%   triangulo cortado introduciria un borde falso justo cerca del haz.
%   Caso de uso principal: eliminar la envolvente lejana (basura del
%   modelo) del cono del reactor, conservando el nucleo del haz.
%
% Ver tambien: loadReactorSTL, signedDistanceMap3D.

if nargin < 2 || isempty(box) || numel(box) ~= 6
    error('BNCT:clipMeshBox:cajaInvalida', ...
        'box debe ser [xmin xmax ymin ymax zmin zmax] (6 valores).');
end
box = double(box(:).');

if isfield(S, 'faces')
    campoCaras = 'faces';
elseif isfield(S, 'triangles')
    campoCaras = 'triangles';
else
    error('BNCT:clipMeshBox:mallaInvalida', ...
        'La malla debe tener campo .faces (o .triangles) y .points.');
end

pts  = double(S.points);
caras = double(S.(campoCaras));
if size(pts, 2) ~= 3 || size(caras, 2) ~= 3
    error('BNCT:clipMeshBox:mallaInvalida', ...
        'points debe ser Nx3 y faces Mx3.');
end

% Vertices dentro de la caja (las 3 coordenadas a la vez)
dentro = pts(:,1) >= box(1) & pts(:,1) <= box(2) & ...
         pts(:,2) >= box(3) & pts(:,2) <= box(4) & ...
         pts(:,3) >= box(5) & pts(:,3) <= box(6);

% Triangulos con los tres vertices dentro
triOk = all(dentro(caras), 2);
kept  = find(triOk);
caras2 = caras(triOk, :);

% Renumerar puntos: quedarse solo con los usados
usados = unique(caras2(:));
if isempty(usados)
    S2 = struct('points', zeros(0,3), campoCaras, zeros(0,3));
    if strcmp(campoCaras, 'triangles')
        S2.faces = S2.triangles;
    end
    return;
end

nuevoIndice = zeros(max(usados), 1);
nuevoIndice(usados) = 1:numel(usados);
caras2 = nuevoIndice(caras2);
pts2 = pts(usados, :);

S2 = struct();
S2.points = pts2;
S2.(campoCaras) = caras2;
if strcmp(campoCaras, 'triangles')
    S2.faces = caras2;
end
end
