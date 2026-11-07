function [collision, dmin] = collisionCheck(patientStruct, wallStruct, bodyStruct, T, safetyDistance)
%COLLISIONCHECK Verifica colision entre la superficie del paciente y el reactor.
%   [collision, dmin] = COLLISIONCHECK(PATIENTSTRUCT, WALLSTRUCT, BODYSTRUCT,
%                                       T, SAFETYDISTANCE)
%
%   PATIENTSTRUCT : struct con .surface.points (Nx3, mm, marco paciente)
%   WALLSTRUCT    : malla del reactor (points/faces); usa el SDF precomputado
%                   (campos .sdf + .gridX/.gridY/.gridZ de
%                   signedDistanceMap3D) si esta disponible -> rapido.
%   BODYSTRUCT    : malla adicional opcional (puede estar vacia).
%   T             : transformacion rigida 4x4 (paciente -> haz), de rigidTransform.
%   SAFETYDISTANCE: distancia de seguridad en mm (OPCIONAL, default 1 mm).
%
%   SALIDA:
%     dmin      : distancia minima firmada del paciente al reactor (mm).
%                 Positivo = fuera del solido; negativo = dentro (penetracion).
%                 Con SDF precomputado es una distancia firmada; con la
%                 verificacion exacta (pointToMeshDistance) es unsigned.
%     collision : true si dmin < safetyDistance.
%
%   ESTRATEGIA:
%     1) Si la malla trae SDF precomputado (recomendado: ver
%        signedDistanceMap3D) se evalua por interpolacion trilineal -> O(N).
%     2) Si no hay SDF y la malla es chica (<5000 triangulos) se usa la
%        distancia EXACTA punto-malla (pointToMeshDistance).
%     3) Si no hay SDF y la malla es grande se usa igualmente la distancia
%        exacta con un aviso (una sola vez por sesion): precalcule el SDF
%        para evaluar miles de poses rapidamente.
%
%   Los puntos que caen fuera del dominio del SDF se evaluan como la
%   distancia a la caja de la grilla SDF (cota inferior conservadora:
%   nunca subestima el riesgo; puede generar colisiones espurias solo en
%   poses extremas). Ver README (limitaciones).

if nargin < 5 || isempty(safetyDistance)
    safetyDistance = 1;   % mm
end

pts = patientStruct.surface.points;
tpts = (T(1:3,1:3) * pts.' + T(1:3,4)).';   % superficie en el marco del haz

d = distanciaReactor(wallStruct, tpts);
dmin = min(d);
if esMallaValida(bodyStruct)
    dmin = min(dmin, min(distanciaReactor(bodyStruct, tpts)));
end

collision = dmin < safetyDistance;
end

% ------------------------------------------------------------------------
function ok = esMallaValida(S)
ok = isstruct(S) && isfield(S, 'points') && isfield(S, 'faces') ...
     && ~isempty(S.faces) && ~isempty(S.points);
end

% ------------------------------------------------------------------------
function d = distanciaReactor(S, pts)
% Distancia de cada punto transformado a la malla del reactor.
if ~esMallaValida(S)
    d = inf;
    return
end

if tieneSDF(S)
    d = evaluaSDF(S, pts);
    return
end

% sin SDF: distancia exacta (aviso una vez por sesion: ~2 ordenes mas lenta)
persistent avisadoSinSDF
if isempty(avisadoSinSDF)
    avisadoSinSDF = false;
end
if ~avisadoSinSDF
    avisadoSinSDF = true;
    warning('BNCT:collisionCheck:sinSDF', ...
        ['Malla de %d triangulos sin SDF precomputado: collisionCheck usa ' ...
         'distancia exacta punto-malla (ordenes de magnitud mas lenta que el ' ...
         'SDF). Precalcule con signedDistanceMap3D y adjunte el resultado a la ' ...
         'malla (campos .sdf + .gridX/.gridY/.gridZ).'], size(S.faces,1));
end
[d, ~] = pointToMeshDistance(pts, S.points, S.faces);
d = min(d);   % distancia exacta es unsigned; sin signo de penetracion
end

% ------------------------------------------------------------------------
function ok = tieneSDF(S)
ok = isfield(S, 'sdf') && isfield(S, 'gridX') && isfield(S, 'gridY') ...
     && isfield(S, 'gridZ') && ~isempty(S.sdf);
end

% ------------------------------------------------------------------------
function d = evaluaSDF(S, pts)
% Evalua el SDF (mm, positivo afuera) por interpolacion trilineal ndgrid.
n = size(pts, 1);
d = inf(n, 1);

inDom = pts(:,1) >= S.gridX(1)  & pts(:,1) <= S.gridX(end) & ...
        pts(:,2) >= S.gridY(1)  & pts(:,2) <= S.gridY(end) & ...
        pts(:,3) >= S.gridZ(1)  & pts(:,3) <= S.gridZ(end);

if any(inDom)
    d(inDom) = interpn(S.gridX, S.gridY, S.gridZ, S.sdf, ...
                       pts(inDom,1), pts(inDom,2), pts(inDom,3), 'linear');
end

% fuera de dominio: distancia a la caja del grid (cota inferior conservadora)
out = ~inDom;
if any(out)
    P  = pts(out, :);
    lo = [S.gridX(1),  S.gridY(1),  S.gridZ(1)];
    hi = [S.gridX(end), S.gridY(end), S.gridZ(end)];
    gap = max(max(lo - P, P - hi), [], 2);
    gap = max(gap, 0);
    d(out) = vecnorm(gap, 2, 2);
end
end
