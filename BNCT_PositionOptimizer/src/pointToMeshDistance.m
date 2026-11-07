function [dmin, idxTri] = pointToMeshDistance(points, meshPoints, meshFaces)
%POINTTOMESHDISTANCE Distancia EXACTA punto-a-malla (punto-triangulo), en mm.
%   [dmin, idxTri] = POINTTOMESHDISTANCE(POINTS, MESHPOINTS, MESHFACES)
%
%   POINTS     : Nx3 coordenadas de consulta (mm).
%   MESHPOINTS : Vx3 vertices de la malla (mm).
%   MESHFACES  : Mx3 indices de triangulo, 1-based sobre MESHPOINTS.
%
%   SALIDA:
%     dmin   : Nx1 distancia minima sin signo punto-superficie (mm).
%     idxTri : Nx1 indice (1-based) del triangulo mas cercano.
%
%   USO: verificacion DETERMINISTA de finalistas (validatePHITS). El SDF
%   voxelizado (signedDistanceMap3D) es rapido pero con precision ~1 voxel;
%   esta funcion es exacta (algoritmo de closest-point-on-triangle de
%   Ericson, "Real-Time Collision Detection"), calculada por chunks para
%   mallas grandes. Coste ~ O(N x M): para mallas de millones de
%   triangulos puede tardar minutos; es un chequeo puntual, no del bucle
%   del optimizador.
%
%   La distancia es SIN signo: no distingue estar dentro/fuera del solido.
%   Para el signo use el SDF; para el chequeo de colision de finalistas
%   combine ambos (dmin exacta < safetyDistance => colision geometrica).

n = size(points, 1);
m = size(meshFaces, 1);

if size(points, 2) ~= 3
    error('BNCT:pointToMeshDistance:puntosMalFormados', ...
        'points debe ser Nx3; recibido %dx%d.', n, size(points, 2));
end
if isempty(meshFaces)
    dmin   = inf(n, 1);
    idxTri = zeros(n, 1);
    return
end
if max(meshFaces(:)) > size(meshPoints,1) || min(meshFaces(:)) < 1
    error('BNCT:pointToMeshDistance:indicesInvalidos', ...
        'meshFaces referencia indices fuera de [1..V] de meshPoints.');
end

points = double(points);
meshPoints = double(meshPoints);
meshFaces = double(meshFaces);

A = meshPoints(meshFaces(:,1), :);           % Mx3
B = meshPoints(meshFaces(:,2), :);
C = meshPoints(meshFaces(:,3), :);
cent = (A + B + C) / 3;
rad = sqrt(max(max(sum((A - cent).^2, 2), sum((B - cent).^2, 2)), ...
                sum((C - cent).^2, 2)));     % radio de esfera envolvente

dmin   = inf(n, 1);
idxTri = zeros(n, 1);

for t = 1:m
    % --- prefiltro por esfera envolvente (cota inferior de la distancia) ---
    dx = points - cent(t, :);
    d2c = sum(dx.^2, 2);
    lb = (sqrt(max(d2c, 0)) - rad(t)).^2;    % (|p-c| - r)^2 <= d(p,tri)^2
    cand = lb < dmin.^2;
    if ~any(cand)
        continue
    end

    P = points(cand, :);                     % solo candidatos
    dist = distPuntoTriangulo(P, A(t,:), B(t,:), C(t,:));

    upd = dist < dmin(cand);
    if any(upd)
        filas = find(cand);
        dmin(filas(upd))   = dist(upd);
        idxTri(filas(upd)) = t;
    end
end
end

% ========================================================================
function dist = distPuntoTriangulo(P, a, b, c)
% Distancia de cada fila de P (Nx3) al triangulo (a,b,c) (1x3 cada uno).
% Implementacion vectorizada de Ericson, "Real-Time Collision Detection"
% (closest point on triangle, regiones de Voronoi 1-7).
ab = b - a;
ac = c - a;

ap = P - a;
d1 = ap * ab.';       % Nx1
d2 = ap * ac.';       % Nx1
bp = P - b;
d3 = bp * ab.';
d4 = bp * ac.';
cp = P - c;
d5 = cp * ab.';
d6 = cp * ac.';

vc = d1 .* d4 - d3 .* d2;
vb = d5 .* d2 - d1 .* d6;
va = d3 .* d6 - d5 .* d4;

n = size(P, 1);
Q = zeros(n, 3);
done = false(n, 1);

% region 1: vertice A
m = ~done & d1 <= 0 & d2 <= 0;
Q(m, :) = repmat(a, nnz(m), 1);
done(m) = true;

% region 2: vertice B
m = ~done & d3 >= 0 & d4 <= d3;
Q(m, :) = repmat(b, nnz(m), 1);
done(m) = true;

% region 3: arista AB (v = d1/(d1-d3))
m = ~done & vc <= 0 & d1 >= 0 & d3 <= 0;
if any(m)
    k = nnz(m);
    den = d1(m) - d3(m);
    v = zeros(k, 1);
    ok = den ~= 0;
    v(ok) = d1(m);
    v(ok) = v(ok) ./ den(ok);
    Q(m, :) = repmat(a, k, 1) + v * ab;
    done(m) = true;
end

% region 4: vertice C
m = ~done & d6 >= 0 & d5 <= d6;
Q(m, :) = repmat(c, nnz(m), 1);
done(m) = true;

% region 5: arista AC (w = d2/(d2-d6))
m = ~done & vb <= 0 & d2 >= 0 & d6 <= 0;
if any(m)
    k = nnz(m);
    den = d2(m) - d6(m);
    w = zeros(k, 1);
    ok = den ~= 0;
    w(ok) = d2(m);
    w(ok) = w(ok) ./ den(ok);
    Q(m, :) = repmat(a, k, 1) + w * ac;
    done(m) = true;
end

% region 6: arista BC (w = (d4-d3)/((d4-d3)+(d5-d6)))
m = ~done & va <= 0 & (d4 - d3) >= 0 & (d5 - d6) >= 0;
if any(m)
    k = nnz(m);
    den = (d4(m) - d3(m)) + (d5(m) - d6(m));
    w = zeros(k, 1);
    ok = den ~= 0;
    w(ok) = d4(m) - d3(m);
    w(ok) = w(ok) ./ den(ok);
    Q(m, :) = repmat(b, k, 1) + w * (c - b);
    done(m) = true;
end

% region 7: interior (coordenadas baricentricas)
m = ~done;
if any(m)
    k = nnz(m);
    den = va(m) + vb(m) + vc(m);
    v = zeros(k, 1);
    w = zeros(k, 1);
    ok = den ~= 0;
    v(ok) = vb(m);
    v(ok) = v(ok) ./ den(ok);
    w(ok) = vc(m);
    w(ok) = w(ok) ./ den(ok);
    Q(m, :) = repmat(a, k, 1) + v * ab + w * ac;
    % degenerados (den==0): vertice mas cercano
    idx = find(m);
    mal = idx(~ok);
    if ~isempty(mal)
        Pm = P(mal, :);
        dA = vecnorm(Pm - a, 2, 2);
        dB = vecnorm(Pm - b, 2, 2);
        dC = vecnorm(Pm - c, 2, 2);
        [~, cual] = min([dA dB dC], [], 2);
        mA = mal(cual == 1);  Q(mA, :) = repmat(a, numel(mA), 1);
        mB = mal(cual == 2);  Q(mB, :) = repmat(b, numel(mB), 1);
        mC = mal(cual == 3);  Q(mC, :) = repmat(c, numel(mC), 1);
    end
end

dist = vecnorm(P - Q, 2, 2);
end
