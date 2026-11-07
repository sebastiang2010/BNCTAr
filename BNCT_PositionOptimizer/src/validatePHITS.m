function [report, bundleDir] = validatePHITS(finalistas, data)
%VALIDATEPHITS Verificacion determinista y robustez de poses finalistas.
%   [report, bundleDir] = VALIDATEPHITS(FINALISTAS, DATA)
%
%   FINALISTAS : matriz Kx6 (una pose por fila, mm/grados), cell array de
%                1x6, o struct array con campo .p.
%   DATA       : mismo struct que optimizePosition (patientStruct,
%                wallStruct, bodyStruct, doseData, radioBio,
%                safetyDistance, options).
%
%   POR CADA FINALISTA:
%     1) Colision EXACTA: pointToMeshDistance de la superficie transformada
%        contra la malla del reactor (independiente del SDF voxelizado).
%     2) Metricas TCP/NTCP/UTCP (via objectiveFunction, SDF rapido).
%     3) Analisis de ROBUSTEZ: perturbaciones +-dx (default 2 mm) por eje
%        de traslacion y +-dtheta (default 2 grados) por eje de rotacion
%        (12 poses + la nominal); se reporta spread de UTCP y cuantas
%        perturbaciones violan limites (NTCP o colision SDF).
%     4) EXPORTACION del bundle para recalculo PHITS (Nivel B):
%        - transform.json  : p, T (4x4), metricas y robustez por finalista
%        - mascaras transformadas en NIfTI (tumor + OARs) via niftiwrite
%          (si no esta disponible, .mat equivalente)
%        - resumen.csv con todas las metricas
%
%   SALIDA:
%     report    : struct array (uno por finalista) con todos los campos.
%     bundleDir : ruta de la carpeta del bundle exportado.
%
%   NIVEL B: el bundle NO ejecuta PHITS; deja todo lo necesario (pose, T,
%   anatomia transformada) para que el usuario lance el recalculo con los
%   modelos reales (cambiar composicion/boro invalida la dosis del Nivel A).

if nargin < 2
    error('BNCT:validatePHITS:faltanArgs', 'Se requieren (finalistas, data).');
end
P = normalizarFinalistas(finalistas);
K = size(P, 1);

opc = struct();
if isfield(data, 'options') && isstruct(data.options)
    opc = data.options;
end
dx     = valor(opc, 'robustDx', 2);        % mm
dtheta = valor(opc, 'robustDtheta', 2);    % grados
bundleDir = valor(opc, 'bundleDir', ...
    fullfile(tempdir, ['BNCT_validacion_', char(datetime('now','Format','yyyyMMdd_HHmmss'))]));
if ~exist(bundleDir, 'dir')
    mkdir(bundleDir);
end

deltas = [dx dx dx dtheta dtheta dtheta];
safety = 1;
if isfield(data, 'safetyDistance') && ~isempty(data.safetyDistance)
    safety = double(data.safetyDistance);
end

report = struct('id', {}, 'p', {}, 'T', {}, 'tcp', {}, 'ntcp', {}, 'utcp', {}, ...
    'dminSDF', {}, 'dminExacta', {}, 'colisionExacta', {}, ...
    'utcpPertMin', {}, 'utcpPertMax', {}, 'utcpPertMedia', {}, ...
    'utcpPertStd', {}, 'nPertViolan', {}, 'archivos', {});

filasCSV = cell(K, 15);
json = struct('generado', char(datetime('now')), 'unidadTraslacion', 'mm', ...
    'unidadRotacion', 'grados', ...
    'convencion', 'T mapea paciente->haz; p=0 es la pose de referencia PHITS', ...
    'finalistas', struct([]));

for k = 1:K
    p = P(k, :);
    T = rigidTransform(p);

    % ---- 1) colision exacta contra las mallas del reactor ----
    ptsHaz = (T(1:3,1:3) * data.patientStruct.surface.points.' + T(1:3,4)).';
    dminExacta = inf;
    if isfield(data, 'wallStruct') && isfield(data.wallStruct, 'faces') ...
            && ~isempty(data.wallStruct.faces)
        [d, ~] = pointToMeshDistance(ptsHaz, data.wallStruct.points, data.wallStruct.faces);
        dminExacta = min(dminExacta, min(d));
    end
    if isfield(data, 'bodyStruct') && isfield(data.bodyStruct, 'faces') ...
            && ~isempty(data.bodyStruct.faces)
        [d, ~] = pointToMeshDistance(ptsHaz, data.bodyStruct.points, data.bodyStruct.faces);
        dminExacta = min(dminExacta, min(d));
    end
    colisionExacta = dminExacta < safety;

    % ---- 2) metricas (mismo camino que el optimizador) ----
    [~, info] = objectiveFunction(p, data);

    % ---- 3) robustez: 12 perturbaciones axiales + nominal ----
    utcpPert = zeros(12, 1);
    violan = 0;
    j = 0;
    for eje = 1:6
        for sgn = [-1, 1]
            j = j + 1;
            q = p;
            q(eje) = q(eje) + sgn * deltas(eje);
            [~, infoPert] = objectiveFunction(q, data);
            utcpPert(j) = infoPert.utcp;
            if ~isfinite(infoPert.utcp) || infoPert.limitsViolated || infoPert.collision
                violan = violan + 1;
            end
        end
    end
    utcpFin = utcpPert(isfinite(utcpPert));
    if isempty(utcpFin)
        utcpFin = -inf;
    end

    % ---- 4) exportar mascaras transformadas ----
    archivos = exportarMascaras(k, p, T, data, bundleDir);

    % ---- acumular reporte ----
    report(k).id             = k;
    report(k).p              = p;
    report(k).T              = T;
    report(k).tcp            = info.tcp;
    report(k).ntcp           = info.ntcp;
    report(k).utcp           = info.utcp;
    report(k).dminSDF        = info.dmin;
    report(k).dminExacta     = dminExacta;
    report(k).colisionExacta = colisionExacta;
    report(k).utcpPertMin    = min(utcpFin);
    report(k).utcpPertMax    = max(utcpFin);
    report(k).utcpPertMedia  = mean(utcpFin);
    report(k).utcpPertStd    = std(utcpFin);
    report(k).nPertViolan    = violan;
    report(k).archivos       = archivos;

    filasCSV(k, :) = {k, p(1), p(2), p(3), p(4), p(5), p(6), ...
        info.tcp, info.ntcp.mucosa, info.ntcp.piel, info.utcp, ...
        dminExacta, min(utcpFin), max(utcpFin), violan};

    fj = struct('id', k, 'p', p, 'T', T, 'tcp', info.tcp, ...
        'ntcpMucosa', info.ntcp.mucosa, 'ntcpPiel', info.ntcp.piel, ...
        'utcp', info.utcp, 'dminSDF_mm', info.dmin, ...
        'dminExacta_mm', dminExacta, 'colisionExacta', colisionExacta, ...
        'utcpPertMin', min(utcpFin), 'utcpPertMax', max(utcpFin), ...
        'utcpPertMedia', mean(utcpFin), 'nPertViolan', violan);
    if isempty(json.finalistas)
        json.finalistas = fj;
    else
        json.finalistas(end+1) = fj;
    end
end

% ---------------- escribir transform.json ----------------
fid = fopen(fullfile(bundleDir, 'transform.json'), 'w');
if fid == -1
    error('BNCT:validatePHITS:noSePuedeEscribir', ...
        'No se pudo crear transform.json en %s.', bundleDir);
end
fwrite(fid, jsonencode(json));
fclose(fid);

% ---------------- escribir resumen.csv ----------------
tbl = array2table(filasCSV, 'VariableNames', {'id', 'x_mm', 'y_mm', 'z_mm', ...
    'alpha_deg', 'beta_deg', 'gamma_deg', 'tcp', 'ntcp_mucosa', 'ntcp_piel', ...
    'utcp', 'dmin_exacta_mm', 'utcp_pert_min', 'utcp_pert_max', 'n_pert_violan'});
writetable(tbl, fullfile(bundleDir, 'resumen.csv'));

for i = 1:numel(report)
    report(i).bundleDir = bundleDir;   % campo de conveniencia en cada finalista
end
end

% ========================================================================
function P = normalizarFinalistas(finalistas)
if isnumeric(finalistas)
    P = double(finalistas);
    if size(P, 2) ~= 6
        error('BNCT:validatePHITS:formatoFinalistas', ...
            'La matriz de finalistas debe ser Kx6.');
    end
    return
end
if iscell(finalistas)
    P = zeros(numel(finalistas), 6);
    for i = 1:numel(finalistas)
        P(i, :) = finalistas{i}(:).';
    end
    return
end
if isstruct(finalistas) && all(isfield(finalistas, 'p'))
    P = zeros(numel(finalistas), 6);
    for i = 1:numel(finalistas)
        P(i, :) = finalistas(i).p(:).';
    end
    return
end
error('BNCT:validatePHITS:formatoFinalistas', ...
    'FINALISTAS debe ser Kx6, cell de 1x6 o struct array con campo .p.');
end

% ========================================================================
function archivos = exportarMascaras(k, p, T, data, bundleDir)
% Exporta tumor + OARs transformados (warp inverso con interpolacion nearest).
nombresOAR = {};
if isfield(data.patientStruct, 'oars')
    nombresOAR = fieldnames(data.patientStruct.oars);
end
archivos = strings(1 + numel(nombresOAR) + 1, 1);  % tumor + OARs + pose.txt
n = 0;
Tinv = eye(4);
Tinv(1:3,1:3) = T(1:3,1:3)';                 % R' (rotacion ortonormal)
Tinv(1:3,4)   = -T(1:3,1:3)' * T(1:3,4);     % -R'*t

% tumor
n = n + 1;
archivos(n) = exportarUna(k, 'tumor', data.patientStruct.tumorMask, ...
    obtenerGridDe(data.patientStruct, 'tumor'), Tinv, bundleDir);

% OARs
for i = 1:numel(nombresOAR)
    entrada = data.patientStruct.oars.(nombresOAR{i});
    if isstruct(entrada) && isfield(entrada, 'mask')
        mascara = entrada.mask;
    else
        mascara = entrada;
    end
    n = n + 1;
    archivos(n) = exportarUna(k, nombresOAR{i}, mascara, ...
        obtenerGridDe(data.patientStruct, nombresOAR{i}), Tinv, bundleDir);
end

% registro de la pose en texto plano (util para 3D Slicer / PHITS)
n = n + 1;
archivos(n) = string(sprintf('finalista%d_pose.txt', k));
fid = fopen(fullfile(bundleDir, archivos(n)), 'w');
fprintf(fid, 'p = [%.6f %.6f %.6f %.6f %.6f %.6f]  %% [x y z] mm, [alpha beta gamma] deg\n', p);
fprintf(fid, 'T = [%.8f %.8f %.8f %.8f;\n', T(1,:));
fprintf(fid, '     %.8f %.8f %.8f %.8f;\n', T(2,:));
fprintf(fid, '     %.8f %.8f %.8f %.8f;\n', T(3,:));
fprintf(fid, '     %.8f %.8f %.8f %.8f];  %% mapea paciente -> haz\n', T(4,:));
fclose(fid);
end

function a = exportarUna(k, nombre, mascara, grids, Tinv, bundleDir)
[X, Y, Z] = ndgrid(grids.gridX, grids.gridY, grids.gridZ);
Pq = [X(:), Y(:), Z(:)];
Pt = Pq * Tinv(1:3,1:3).' + Tinv(1:3,4).';   % expansion implicita (Nx3 + 1x3)
val = interpn(grids.gridX, grids.gridY, grids.gridZ, double(mascara > 0), ...
    Pt(:,1), Pt(:,2), Pt(:,3), 'nearest');
val(isnan(val)) = 0;   % fuera del dominio de la grilla -> fondo
maskT = double(reshape(val > 0.5, size(mascara)));

base = sprintf('finalista%d_%s', k, nombre);
try
    niftiwrite(maskT, fullfile(bundleDir, [base, '.nii']));   % IPT
    a = string([base, '.nii']);
catch
    rutaMat = fullfile(bundleDir, [base, '.mat']);
    save(rutaMat, 'maskT');
    a = string([base, '.mat']);
end
end

function g = obtenerGridDe(ps, nombreRegion)
g = struct('gridX', ps.gridX, 'gridY', ps.gridY, 'gridZ', ps.gridZ);
if strcmp(nombreRegion, 'tumor') && isfield(ps, 'tumorGrid') && isstruct(ps.tumorGrid)
    if isfield(ps.tumorGrid, 'gridX'); g.gridX = ps.tumorGrid.gridX; end
    if isfield(ps.tumorGrid, 'gridY'); g.gridY = ps.tumorGrid.gridY; end
    if isfield(ps.tumorGrid, 'gridZ'); g.gridZ = ps.tumorGrid.gridZ; end
elseif ~strcmp(nombreRegion, 'tumor') && isfield(ps, 'oars') ...
        && isfield(ps.oars, nombreRegion) && isstruct(ps.oars.(nombreRegion))
    src = ps.oars.(nombreRegion);
    if isfield(src, 'gridX'); g.gridX = src.gridX; end
    if isfield(src, 'gridY'); g.gridY = src.gridY; end
    if isfield(src, 'gridZ'); g.gridZ = src.gridZ; end
end
end

function v = valor(S, campo, def)
if isfield(S, campo) && ~isempty(S.(campo))
    v = S.(campo);
else
    v = def;
end
end
