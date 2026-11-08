function sdfStruct = signedDistanceMap3D(meshStruct, options)
%SIGNEDDISTANCEMAP3D Campo de distancia con signo (SDF) del reactor en grilla de voxels.
%   sdfStruct = SIGNEDDISTANCEMAP3D(MESHSTRUCT, OPTIONS)
%
%   MESHSTRUCT : struct con .points (Nx3, mm) y .faces (Mx3, 1-based).
%   OPTIONS    : struct opcional con:
%     .resolution : resolucion isotropa de la grilla en mm (default 2).
%     .margin     : margen alrededor del bbox de la malla, mm. Escalar o
%                   1x3 por eje (default 50). DEBE cubrir toda la anatomia
%                   del paciente y su rango de movimiento.
%     .gridX/.gridY/.gridZ : ejes ndgrid (mm) explicitos (overrides).
%     .fillMethod : 'parity' (default) relleno solido por test de paridad
%                   por columna z; 'shell' trata la superficie rasterizada
%                   (dilatada 1 voxel) como solido de pared fina.
%
%   SALIDA: struct con
%     .sdf       : array 3D [nx ny nz] en mm, POSITIVO AFUERA del solido,
%                  negativo adentro (sdf = bwdist(occ) - bwdist(~occ)).
%     .gridX/.gridY/.gridZ : ejes ndgrid en mm (dim1=x, dim2=y, dim3=z).
%     .voxelSize : [1x3] mm (isotropo = resolution).
%     .points/.faces      : copia de la malla de entrada.
%     .occupancy : logico 3D del solido rasterizado (debug).
%     .fillMethod, .resolution
%
%   POR QUE NO alphaShape: alphaShape cierra concavidades (la abertura del
%   puerto quedaria tapada). Aqui el solido se reconstruye por paridad de
%   cruces por columna (ray casting a lo largo de z), que PRESERVA tuneles
%   perpendiculares a z (la abertura del haz). Limitaciones (README):
%     - Tuneles PARALELOS al eje de relleno (z) se cerrarian de forma
%       erronea: use .fillMethod='shell' o reoriente el relleno.
%     - Precision ~1 voxel (resolution mm). El SDF es para cribado rapido
%       durante la optimizacion; la verificacion determinista de finalistas
%       se hace con pointToMeshDistance (exacta).
%     - 'parity' requiere malla cerrada (2-manifold).
%
%   Requiere Image Processing Toolbox (bwdist; imdilate solo en 'shell').

if nargin < 2 || isempty(options)
    options = struct();
end

% ------------------------- opciones -------------------------
res = 2;   % mm
if isfield(options, 'resolution') && ~isempty(options.resolution)
    res = double(options.resolution);
end
if res <= 0
    error('BNCT:signedDistanceMap3D:resolucionInvalida', ...
        'options.resolution debe ser > 0 mm.');
end

fillMethod = 'parity';
if isfield(options, 'fillMethod') && ~isempty(options.fillMethod)
    fillMethod = lower(options.fillMethod);
end
if ~ismember(fillMethod, {'parity', 'shell'})
    error('BNCT:signedDistanceMap3D:fillMethodInvalido', ...
        'options.fillMethod debe ser ''parity'' o ''shell'' (recibido %s).', fillMethod);
end

pts = double(meshStruct.points);
fcs = double(meshStruct.faces);
if isempty(fcs)
    error('BNCT:signedDistanceMap3D:mallaVacia', ...
        'La malla del reactor no tiene triangulos.');
end
if max(fcs(:)) > size(pts,1) || min(fcs(:)) < 1
    error('BNCT:signedDistanceMap3D:indicesInvalidos', ...
        'faces referencia indices fuera de rango de points.');
end

% ------------------------- grilla -------------------------
if all(isfield(options, {'gridX','gridY','gridZ'})) && ...
        ~isempty(options.gridX) && ~isempty(options.gridY) && ~isempty(options.gridZ)
    gridX = double(options.gridX(:)).';
    gridY = double(options.gridY(:)).';
    gridZ = double(options.gridZ(:)).';
else
    margen = [50 50 50];   % mm
    if isfield(options, 'margin') && ~isempty(options.margin)
        m = double(options.margin(:)).';
        if isscalar(m)
            margen = m * [1 1 1];
        elseif numel(m) == 3
            margen = m;
        end
    end
    lo = min(pts, [], 1) - margen;
    hi = max(pts, [], 1) + margen;
    n = max(2, ceil((hi - lo) / res) + 1);
    gridX = lo(1) + (0:n(1)-1) * res;
    gridY = lo(2) + (0:n(2)-1) * res;
    gridZ = lo(3) + (0:n(3)-1) * res;
end
nx = numel(gridX); ny = numel(gridY); nz = numel(gridZ);
origin = [gridX(1), gridY(1), gridZ(1)];

% ------------------------- rasterizar superficie -------------------------
switch fillMethod
    case 'parity'
        % Malla cerrada: rasterizado por columnas + relleno solido por paridad
        occ = rasterizarSuperficie(pts, fcs, origin, res, [nx ny nz]);
        occ = rellenarParidad(pts, fcs, gridX, gridY, gridZ, origin, res, occ);
    case 'shell'
        % Parche ABIERTO (p. ej. cono del haz sin tapa en la boca x=0): no
        % hay clasificacion de interior. Rasterizado FINO: solo voxels cuyo
        % centro dista <= res/2 del plano de algun triangulo, para no
        % engordar la banda cuando los triangulos son grandes u oblicuos.
        occ = rasterizarSuperficieFina(pts, fcs, origin, res, [nx ny nz]);
        occ = imdilate(occ, ones(3,3,3));   % pared de ~1-2 voxels de espesor
end

% ------------------------- SDF (positivo afuera) -------------------------
% bwdist devuelve distancias en unidades de voxel; con grilla isotropa
% basta multiplicar por res para obtener milimetros.
sdf = res * (bwdist(occ) - bwdist(~occ));

sdfStruct = struct();
sdfStruct.sdf        = sdf;
sdfStruct.gridX      = gridX;
sdfStruct.gridY      = gridY;
sdfStruct.gridZ      = gridZ;
sdfStruct.voxelSize  = [res res res];
sdfStruct.points     = pts;
sdfStruct.faces      = fcs;
sdfStruct.occupancy  = occ;
sdfStruct.fillMethod = fillMethod;
sdfStruct.resolution = res;
end

% ========================================================================
function occ = rasterizarSuperficie(pts, fcs, origin, res, sizeGrid)
% Marca true en los voxels cuyo centro cae dentro de algun triangulo.
% Para cada triangulo se proyecta al plano perpendicular al eje dominante
% de su normal y se hace el test baricentrico 2D sobre los centros de
% voxel del bbox; a lo largo del eje dominante se marca todo el rango.
nx = sizeGrid(1); ny = sizeGrid(2); nz = sizeGrid(3);
occ = false(nx, ny, nz);
escala = norm(max(pts, [], 1) - min(pts, [], 1));   % mm
umbral = 1e-12 * escala^2;

nT = size(fcs, 1);
for t = 1:nT
    v1 = pts(fcs(t,1), :);
    v2 = pts(fcs(t,2), :);
    v3 = pts(fcs(t,3), :);

    ijk  = ([v1; v2; v3] - origin) / res + 1;      % bbox en indices
    iMin = max(1, floor(min(ijk, [], 1)));
    iMax = min(sizeGrid(:).', ceil(max(ijk, [], 1)));
    if any(iMin > iMax)
        continue
    end

    nrm = cross(v2 - v1, v3 - v1);
    [~, ax] = max(abs(nrm));
    uv = setdiff(1:3, ax);                         % ejes de proyeccion

    [Iu, Iv] = ndgrid(iMin(uv(1)):iMax(uv(1)), iMin(uv(2)):iMax(uv(2)));
    uq = origin(uv(1)) + (Iu - 1) * res;           % centros de voxel (mm)
    vq = origin(uv(2)) + (Iv - 1) * res;

    u1 = v1(uv(1)); v1c = v1(uv(2));
    u2 = v2(uv(1)); v2c = v2(uv(2));
    u3 = v3(uv(1)); v3c = v3(uv(2));

    det = (u2 - u1) * (v3c - v1c) - (u3 - u1) * (v2c - v1c);
    dentro = false(size(Iu));
    if abs(det) > umbral
        e1 = (u2 - u1) .* (vq - v1c) - (v2c - v1c) .* (uq - u1);
        e2 = (u3 - u2) .* (vq - v2c) - (v3c - v2c) .* (uq - u2);
        e3 = (u1 - u3) .* (vq - v3c) - (v1c - v3c) .* (uq - u3);
        s = sign(det);
        dentro = (s*e1 >= 0) & (s*e2 >= 0) & (s*e3 >= 0);
    end

    if any(dentro(:))
        IuM = Iu(dentro);
        IvM = Iv(dentro);
        rangoDom = iMin(ax):iMax(ax);
        nxny = nx * ny;
        for w = rangoDom
            switch ax
                case 1    % eje dominante x: (Iu,Iv) = (j,k)
                    lin = w + (IuM - 1) * nx + (IvM - 1) * nxny;
                case 2    % dominante y: (Iu,Iv) = (i,k)
                    lin = IuM + (w - 1) * nx + (IvM - 1) * nxny;
                otherwise % dominante z: (Iu,Iv) = (i,j)
                    lin = IuM + (IvM - 1) * nx + (w - 1) * nxny;
            end
            occ(lin) = true;
        end
    end

    % garantizar presencia de los vertices (triangulos finos/degenerados)
    vertices = [v1; v2; v3];
    for q = 1:3
        idx = min(max(round((vertices(q,:) - origin) / res + 1), 1), sizeGrid(:).');
        occ(idx(1), idx(2), idx(3)) = true;
    end
end
end

% ========================================================================
function occ = rellenarParidad(pts, fcs, gridX, gridY, gridZ, origin, res, occ)
% Solidifica por paridad de cruces por columna (i,j) a lo largo de z:
% para cada columna se ordenan las alturas de cruce con los triangulos y
% se rellena entre pares (dentro-fuera-dentro...). Preserva tuneles con
% eje NO paralelo a z (como la abertura del haz a lo largo de x).
nx = numel(gridX); ny = numel(gridY);
escala = norm(max(pts, [], 1) - min(pts, [], 1));
umbral = 1e-12 * escala^2;

hitsIc = cell(1, size(fcs,1));
hitsJc = cell(1, size(fcs,1));
hitsZc = cell(1, size(fcs,1));

for t = 1:size(fcs, 1)
    v1 = pts(fcs(t,1), :);
    v2 = pts(fcs(t,2), :);
    v3 = pts(fcs(t,3), :);

    nrm = cross(v2 - v1, v3 - v1);
    if abs(nrm(3)) <= umbral
        continue   % plano (casi) vertical: no aporta cruces en z
    end

    ijk  = ([v1; v2; v3] - origin) / res + 1;
    iMin = max(1, floor(min(ijk, [], 1)));
    iMax = min([nx ny Inf], ceil(max(ijk, [], 1)));
    if iMin(1) > iMax(1) || iMin(2) > iMax(2)
        continue
    end

    [II, JJ] = ndgrid(iMin(1):iMax(1), iMin(2):iMax(2));
    xq = gridX(II);
    yq = gridY(JJ);

    det = (v2(1)-v1(1)) * (v3(2)-v1(2)) - (v3(1)-v1(1)) * (v2(2)-v1(2));
    if abs(det) <= umbral
        continue
    end
    e1 = (v2(1)-v1(1)) .* (yq - v1(2)) - (v2(2)-v1(2)) .* (xq - v1(1));
    e2 = (v3(1)-v2(1)) .* (yq - v2(2)) - (v3(2)-v2(2)) .* (xq - v2(1));
    e3 = (v1(1)-v3(1)) .* (yq - v3(2)) - (v1(2)-v3(2)) .* (xq - v3(1));
    s = sign(det);
    dentro = (s*e1 >= 0) & (s*e2 >= 0) & (s*e3 >= 0);
    if ~any(dentro(:))
        continue
    end

    zh = v1(3) - (nrm(1) * (xq - v1(1)) + nrm(2) * (yq - v1(2))) / nrm(3);
    hitsIc{t} = II(dentro);
    hitsJc{t} = JJ(dentro);
    hitsZc{t} = zh(dentro);
end

hitsI = vertcat(hitsIc{:});
hitsJ = vertcat(hitsJc{:});
hitsZ = vertcat(hitsZc{:});
if isempty(hitsI)
    return
end

% agrupar por columna y rellenar entre pares de cruces
[cols, ~, g] = unique([hitsI, hitsJ], 'rows');
for c = 1:size(cols, 1)
    zs = sort(hitsZ(g == c));
    % deduplicar cruces casi identicos (aristas compartidas / grazing)
    tol = 1e-6 * max(1, max(abs(zs))) + 1e-9;
    keep = [true; diff(zs) > tol];
    zs = zs(keep);
    if numel(zs) < 2
        continue
    end
    if mod(numel(zs), 2) == 1
        zs(end) = [];   % numero impar de cruces: descartar el ultimo
    end
    ii = cols(c, 1);
    jj = cols(c, 2);
    for pr = 1:2:numel(zs)
        k1 = find(gridZ >= zs(pr), 1);
        k2 = find(gridZ <= zs(pr+1), 1, 'last');
        if ~isempty(k1) && ~isempty(k2) && k1 <= k2
            occ(ii, jj, k1:k2) = true;
        end
    end
end
end

function occ = rasterizarSuperficieFina(pts, fcs, origin, res, sizeGrid)
% Rasterizado fino para parches abiertos (fillMethod='shell').
% Marca true SOLO en los voxels cuyo centro dista <= res/2 del plano de
% algun triangulo y cuya proyeccion cae dentro del triangulo (test
% baricentrico 2D sobre el plano perpendicular al eje dominante de la
% normal). A diferencia de rasterizarSuperficie, NO marca el rango
% completo del bbox a lo largo del eje dominante: con triangulos grandes
% u oblicuos (cono real del haz) eso engordaba la banda ocupada decenas
% de voxels y generaba distancias negativas fantasma dentro del tunel.
nx = sizeGrid(1); ny = sizeGrid(2); nz = sizeGrid(3);
occ = false(nx, ny, nz);
escala = norm(max(pts, [], 1) - min(pts, [], 1));   % mm
umbral = 1e-12 * escala^2;
semiVoxel = res / 2;

nT = size(fcs, 1);
for t = 1:nT
    v1 = pts(fcs(t,1), :);
    v2 = pts(fcs(t,2), :);
    v3 = pts(fcs(t,3), :);

    ijk  = ([v1; v2; v3] - origin) / res + 1;      % bbox en indices
    iMin = max(1, floor(min(ijk, [], 1)) - 1);
    iMax = min(sizeGrid(:).', ceil(max(ijk, [], 1)) + 1);
    if any(iMin > iMax)
        continue
    end

    nrm = cross(v2 - v1, v3 - v1);
    n2 = sum(nrm.^2);
    if n2 <= umbral
        continue;                                  % triangulo degenerado
    end
    [~, ax] = max(abs(nrm));
    uv = setdiff(1:3, ax);

    % candidatos: centros de voxel del bbox (con 1 voxel de holgura)
    [I, J, K] = ndgrid(iMin(1):iMax(1), iMin(2):iMax(2), iMin(3):iMax(3));
    centros = origin + ([I(:), J(:), K(:)] - 1) * res;   % mm

    % distancia con signo al plano del triangulo
    dPlano = ((centros - v1) * nrm.') / sqrt(n2);
    cerca = abs(dPlano) <= semiVoxel;
    if ~any(cerca)
        continue
    end

    % test baricentrico 2D (mismo convenio que rasterizarSuperficie)
    uq = centros(cerca, uv(1)); vq = centros(cerca, uv(2));
    u1 = v1(uv(1)); v1c = v1(uv(2));
    u2 = v2(uv(1)); v2c = v2(uv(2));
    u3 = v3(uv(1)); v3c = v3(uv(2));
    det = (u2 - u1) * (v3c - v1c) - (u3 - u1) * (v2c - v1c);
    if abs(det) <= umbral
        continue
    end
    e1 = (u2 - u1) .* (vq - v1c) - (v2c - v1c) .* (uq - u1);
    e2 = (u3 - u2) .* (vq - v2c) - (v3c - v2c) .* (uq - u2);
    e3 = (u1 - u3) .* (vq - v3c) - (v1c - v3c) .* (uq - u3);
    s = sign(det);
    dentro = (s*e1 >= 0) & (s*e2 >= 0) & (s*e3 >= 0);
    if ~any(dentro)
        continue
    end

    lin = sub2ind([nx ny nz], I(cerca), J(cerca), K(cerca));
    occ(lin(dentro)) = true;
end

% garantizar presencia de los vertices (triangulos finos/degenerados)
vertices = pts(unique(fcs(:)), :);
for q = 1:size(vertices, 1)
    idx = min(max(round((vertices(q,:) - origin) / res + 1), 1), sizeGrid(:).');
    occ(idx(1), idx(2), idx(3)) = true;
end
end
