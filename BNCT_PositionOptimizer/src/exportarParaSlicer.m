function [archivos, info] = exportarParaSlicer(data, poses, outDir, options)
%EXPORTARPARASLICER Exporta un bundle visualizable en 3D Slicer (marco RAS).
%   [ARCHIVOS, INFO] = EXPORTARPARASLICER(DATA, POSES, OUTDIR, OPTIONS)
%
%   Genera en OUTDIR un bundle para inspeccion visual en 3D Slicer:
%
%     puerto_reactor.vtk      malla del puerto (DATA.wallStruct), marco RAS
%     paciente_pose_<i>.vtk   paciente con rigidTransform(POSES(i,:)), RAS
%     dosis_<comp>.nii        dosis compuesta en la grilla PHITS (affine RAS)
%     transform.json          convencion RAS, grilla de dosis y poses
%     pose_<i>.txt            pose legible: p, dmin y matriz T 4x4
%
%   CONVENCION RAS (requisito central): 3D Slicer trabaja en RAS, mientras
%   que los VTK y la grilla de dosis del proyecto estan en el MARCO DEL HAZ
%   (igual que los .vtk de entrada de data/). Conversion elegida LPS->RAS:
%   (x,y,z) -> (-x,-y,z).
%     1) OPTIONS.espejarRAS (default TRUE): niega x e y de TODOS los puntos
%        exportados en VTK (puerto y pacientes). El espejo se aplica SOLO a
%        los puntos exportados (DESPUES de la transformacion de pose);
%        nunca modifica DATA de entrada.
%     2) La dosis NIfTI no se remuestrea: el espejo va en el affine sform,
%        que Slicer interpreta como voxel 0-based -> RAS:
%           M = [-sx 0 0 -gridX(1); 0 -sy 0 -gridY(1); 0 0 sz gridZ(1); 0 0 0 1]
%        con sx = espaciado de gridX (promedio si la grilla no es uniforme,
%        con aviso). Asi el voxel (i0,j0,k0) 0-based cae en
%        RAS = (-gridX(i0+1), -gridY(j0+1), gridZ(k0+1)) mm.
%        ORDEN DE EJES: los arrays de dosis son [nx ny nz] en ndgrid
%        (dim1<->x, dim2<->y, dim3<->z; ver loadPHITSDose). NIfTI-1 guarda
%        dim1 como eje de indice mas rapido y el affine mapea dim1->x: el
%        orden ya coincide con NIfTI. NO permutar el array.
%        Verificado en R2023b: niftiwrite no acepta una matriz 4x4 suelta;
%        se escribe el volumen, se lee niftiinfo y se reescribe con
%        info.TransformName='Sform' e info.Transform=affine3d(M.') (nifti-
%        write guarda las FILAS de M en srow_x/y/z y sform_code=1, que es
%        lo que Slicer/ITK leen con precedencia sobre qform).
%     3) transform.json documenta SIEMPRE M_RAS = diag([-1 -1 1 1]) y la
%        nota para regenerar sin espejo (espejarRAS=false) si la escena
%        original de Slicer cargaba los .vtk de entrada sin espejar.
%
%   ENTRADAS:
%     DATA    : struct como el de run_optimizer (makeDemoData o modo real):
%       .patientStruct.surface.points : Nx3 mm, marco paciente [obligatorio]
%       .patientStruct.surface.faces  : Mx3 (opcional; si falta, la
%                                       superficie se exporta como puntos)
%       .wallStruct.points/.faces     : malla del puerto [obligatorio]
%       .wallStruct.sdf + .gridX/.gridY/.gridZ : SDF precomputado (opcional)
%       .bodyStruct                   : malla extra de colision (opcional)
%       .doseData                     : canonico de loadPHITSDose (opcional)
%       .radioBio.rbeTumor/rbeNormal  : pesos [B H N gamma] (opcionales)
%       .safetyDistance               : mm (default 1)
%     POSES   : Kx6 [x y z alfa beta gamma] (mm, grados), una pose por fila;
%               [] usa [0 0 0 0 0 0] (paciente sin transformar). Acepta
%               ademas cell de 1x6 y struct array con campo .p.
%     OUTDIR  : carpeta destino (se crea si no existe; solo se sobreescriben
%               los archivos con el mismo nombre, el resto no se toca).
%               Default: <raiz del proyecto>/slicer_out
%     OPTIONS : struct opcional:
%       .espejarRAS      : true (default) | false
%       .doseComponentes : {'total'} (default); 'tumor'/'normal' ponderan
%                          con radioBio.rbeTumor/rbeNormal (Gy-Eq x
%                          doseScale). Si faltan los pesos se exporta solo
%                          'total' (con aviso).
%       .generarDosis    : true (default) | false
%
%   SALIDAS:
%     ARCHIVOS : cellstr con las rutas completas de los archivos escritos.
%     INFO     : struct con .outDir .espejarRAS .safetyDistance .nPoses
%                .poses (struct array) .doses (struct).
%
%   NOTAS:
%     - dmin por pose: collisionCheck con la data tal cual (distancia FIRMA
%       si hay SDF: negativa = penetracion). Si wallStruct no tiene SDF se
%       usa la distancia EXACTA pointToMeshDistance contra wallStruct y
%       bodyStruct (minimo de ambos, SIN signo), con aviso.
%     - 'total' replica la dosis fisica compuesta que usa el costo: con
%       pesos unitarios, calculateUTCP (llamado por objectiveFunction)
%       compone Deq = [D_B D_H D_N D_gamma] * pesos, es decir la suma simple
%       de las 4 componentes PHITS; se multiplica por doseScale (el loader
%       no lo aplica). 'tumor'/'normal' usan los pesos Gy->Gy-Eq de radioBio
%       (misma composicion lineal que calculateUTCP).
%     - Unidades: mm, grados, Gy. Sin figuras, sin red, sin secretos.
%
%   Ver tambien: rigidTransform, collisionCheck, loadPHITSDose,
%   objectiveFunction, calculateUTCP, validatePHITS, signedDistanceMap3D.

% =========================================================================
% 1) defaults y opciones
% =========================================================================
if nargin < 2 || isempty(poses)
    poses = zeros(1, 6);                       % paciente sin transformar
end
if nargin < 3 || isempty(outDir)
    outDir = fullfile(fileparts(mfilename('fullpath')), '..', 'slicer_out');
end
if nargin < 4 || isempty(options)
    options = struct();
end
if ~isstruct(options)
    error('BNCT:exportarParaSlicer:opcionesInvalidas', ...
        'OPTIONS debe ser un struct (o omitirse).');
end

espejarRAS   = logical(valorConDefault(options, 'espejarRAS', true));
generarDosis = logical(valorConDefault(options, 'generarDosis', true));
componentes  = normalizarComponentes( ...
    valorConDefault(options, 'doseComponentes', {'total'}));

if espejarRAS
    factorEspejo = [-1 -1 1];   % (x,y,z) -> (-x,-y,z): marco del haz -> RAS
else
    factorEspejo = [1 1 1];     % sin espejo: queda el marco del haz
end

% =========================================================================
% 2) validacion de DATA y de las mallas
% =========================================================================
if ~isstruct(data) || ~isfield(data, 'patientStruct') || ...
        ~isfield(data.patientStruct, 'surface') || ...
        ~isfield(data.patientStruct.surface, 'points') || ...
        isempty(data.patientStruct.surface.points)
    error('BNCT:exportarParaSlicer:sinSuperficiePaciente', ...
        ['DATA.patientStruct.surface.points no existe (Nx3, mm, marco ' ...
         'paciente). En modo real se llena desde BODY.vtk; en demo lo arma ' ...
         'makeDemoData (ver run_optimizer.m).']);
end
if ~isfield(data, 'wallStruct') || ~isfield(data.wallStruct, 'points') || ...
        isempty(data.wallStruct.points) || ...
        ~isfield(data.wallStruct, 'faces') || isempty(data.wallStruct.faces)
    error('BNCT:exportarParaSlicer:sinMallaReactor', ...
        ['DATA.wallStruct.points/.faces no existen (malla del puerto, mm, ' ...
         'marco del haz). Use loadReactorSTL + clipMeshBox (ver ' ...
         'run_optimizer.m, modo real).']);
end

ptsPaciente = double(data.patientStruct.surface.points);
if size(ptsPaciente, 2) ~= 3
    error('BNCT:exportarParaSlicer:puntosPacienteMalFormados', ...
        'surface.points debe ser Nx3; recibio %dx%d.', ...
        size(ptsPaciente, 1), size(ptsPaciente, 2));
end
ptsPuerto   = double(data.wallStruct.points);
carasPuerto = double(data.wallStruct.faces);
validarCaras(carasPuerto, size(ptsPuerto, 1), 'wallStruct');

carasPaciente = [];
if isfield(data.patientStruct.surface, 'faces') && ...
        ~isempty(data.patientStruct.surface.faces)
    carasPaciente = double(data.patientStruct.surface.faces);
    validarCaras(carasPaciente, size(ptsPaciente, 1), 'patientStruct.surface');
else
    warning('BNCT:exportarParaSlicer:sinCarasPaciente', ...
        ['patientStruct.surface no trae .faces: paciente_pose_<i>.vtk se ' ...
         'exporta como nube de puntos (VERTICES). Agregue surface.faces ' ...
         '(p.ej. de BODY.vtk) para obtener una superficie cerrada.']);
end

wallStruct = data.wallStruct;
bodyStruct = struct('points', zeros(0, 3), 'faces', zeros(0, 3));
if isfield(data, 'bodyStruct') && isstruct(data.bodyStruct)
    bodyStruct = data.bodyStruct;
end
safetyDistance = 1;                            % mm (default de collisionCheck)
if isfield(data, 'safetyDistance') && ~isempty(data.safetyDistance)
    safetyDistance = double(data.safetyDistance);
end

radioBio = struct();
if isfield(data, 'radioBio') && isstruct(data.radioBio)
    radioBio = data.radioBio;
end
doseData = [];
if isfield(data, 'doseData') && isstruct(data.doseData)
    doseData = data.doseData;
end

P = normalizarPoses(poses);
K = size(P, 1);

% SDF precomputado del puerto: decide el metodo de dmin (una sola vez)
haySDF = tieneSDF(wallStruct);
if ~haySDF
    warning('BNCT:exportarParaSlicer:sinSDF', ...
        ['wallStruct no tiene SDF precomputado: dmin se calcula con ' ...
         'pointToMeshDistance (exacta, SIN signo) contra wallStruct y ' ...
         'bodyStruct. Para distancia firmada (penetracion) precalcule con ' ...
         'signedDistanceMap3D (ver run_optimizer.m).']);
end

% =========================================================================
% 3) carpeta de salida
% =========================================================================
if ~exist(outDir, 'dir')
    [okDir, msgDir] = mkdir(outDir);
    if ~okDir
        error('BNCT:exportarParaSlicer:noSePuedeCrearDir', ...
            'No se pudo crear la carpeta %s: %s', outDir, msgDir);
    end
end

% =========================================================================
% 4) plan de dosis: componentes, grilla y affine RAS
% =========================================================================
[listaDosis, pesosRBE] = planificarDosis(componentes, doseData, radioBio, generarDosis);

M = [];                % affine 4x4 voxel(0-based) -> RAS (pre-multiplicacion)
espaciado = [];
if ~isempty(listaDosis)
    gridX = double(doseData.gridX(:).');
    gridY = double(doseData.gridY(:).');
    gridZ = double(doseData.gridZ(:).');
    espaciado = espaciadoDeGrilla(gridX, gridY, gridZ);   % avisa si no uniforme

    % AFFINE RAS: Slicer lee sform como mapeo de voxel 0-based -> mundo.
    % Con el espejo, RAS_x(i0) = -gridX(1) - i0*sx = -gridX(i0+1) (el
    % indice MATLAB de la dimension 1 es i0+1). Ver comentario del header.
    if espejarRAS
        M = [-espaciado(1), 0, 0, -gridX(1);
             0, -espaciado(2), 0, -gridY(1);
             0, 0, espaciado(3), gridZ(1);
             0, 0, 0, 1];
    else
        M = [espaciado(1), 0, 0, gridX(1);
             0, espaciado(2), 0, gridY(1);
             0, 0, espaciado(3), gridZ(1);
             0, 0, 0, 1];
    end
end

% =========================================================================
% 5) puerto_reactor.vtk (espejo RAS aplicado solo a los puntos exportados)
% =========================================================================
% preallocacion: puerto + 2 archivos por pose + un .nii por componente
archivos = cell(1, 1 + 2 * K + numel(listaDosis));
nArch = 0;
nArch = nArch + 1;
archivos{nArch} = escribirVTK(fullfile(outDir, 'puerto_reactor.vtk'), ...
    ptsPuerto .* factorEspejo, carasPuerto, 'Puerto reactor BNCT (RAS)');

% =========================================================================
% 6) paciente_pose_<i>.vtk + pose_<i>.txt + dmin por pose
% =========================================================================
poseInfo = struct('indice', {}, 'p', {}, 'T', {}, 'dmin_mm', {}, ...
    'colisiona', {}, 'dminMetodo', {}, 'archivoVTK', {}, 'archivoTxt', {});

for i = 1:K
    p = P(i, :);
    T = rigidTransform(p);                     % paciente -> haz (SIN espejo)

    % superficie transformada al marco del haz (mm)
    ptsHaz = (T(1:3, 1:3) * ptsPaciente.' + T(1:3, 4)).';

    % --- dmin: collisionCheck con SDF, o distancia exacta sin SDF ---
    if haySDF
        [colisiona, dmin] = collisionCheck(data.patientStruct, wallStruct, ...
            bodyStruct, T, safetyDistance);
        metodo = 'sdf_firmada';                % negativa = penetracion
    else
        dmin = distanciaExacta(ptsHaz, wallStruct, bodyStruct);
        colisiona = dmin < safetyDistance;
        metodo = 'exacta_sin_signo';
    end

    % --- VTK del paciente (espejo RAS solo en los puntos exportados) ---
    nombreVTK = sprintf('paciente_pose_%d.vtk', i);
    rutaVTK = fullfile(outDir, nombreVTK);
    nArch = nArch + 1;
    archivos{nArch} = escribirVTK(rutaVTK, ptsHaz .* factorEspejo, ...
        carasPaciente, sprintf('Paciente BNCT pose %d (RAS)', i));

    % --- pose_<i>.txt legible ---
    nombreTxt = sprintf('pose_%d.txt', i);
    rutaTxt = fullfile(outDir, nombreTxt);
    escribirPoseTxt(rutaTxt, i, K, p, T, dmin, colisiona, metodo, ...
        safetyDistance);
    nArch = nArch + 1;
    archivos{nArch} = rutaTxt;

    poseInfo(i) = struct('indice', i, 'p', p, 'T', T, 'dmin_mm', dmin, ...
        'colisiona', colisiona, 'dminMetodo', metodo, ...
        'archivoVTK', rutaVTK, 'archivoTxt', rutaTxt);
end

% =========================================================================
% 7) dosis_<comp>.nii (affine RAS; NO se permuta el array, ver header)
% =========================================================================
dosesInfo = struct('exportado', false, 'componentes', {{}}, 'doseScale', 1, ...
    'unidades', '', 'origen_RAS', [], 'spacing_RAS', [], 'dims', [], ...
    'conAffine', false, 'archivos', {{}});
nombresNii = cell(1, numel(listaDosis));

if ~isempty(listaDosis)
    validarDoseData(doseData);                 % grillas y arrays consistentes
    escalaDosis = valorConDefault(doseData, 'doseScale', 1);
    if isempty(escalaDosis) || ~isnumeric(escalaDosis)
        escalaDosis = 1;
    end
    escalaDosis = double(escalaDosis(1));

    for c = 1:numel(listaDosis)
        [D, tipoDosis] = construirDosisCompuesta(listaDosis{c}, doseData, ...
            pesosRBE, escalaDosis);
        nombreNii = sprintf('dosis_%s.nii', lower(listaDosis{c}));
        rutaNii = fullfile(outDir, nombreNii);
        conAffine = escribirNIfTI(D, rutaNii, M, espaciado);
        nArch = nArch + 1;
        archivos{nArch} = rutaNii;
        nombresNii{c} = nombreNii;
        dosesInfo.conAffine = dosesInfo.conAffine || conAffine;
        if c == 1
            dosesInfo.dims = size(D);
            dosesInfo.origen_RAS = M(1:3, 4).';
            dosesInfo.spacing_RAS = espaciado;
        end
        if strcmpi(listaDosis{c}, 'total') || strcmpi(tipoDosis, 'Gy')
            dosesInfo.unidades = 'Gy (doseScale aplicado)';
        else
            dosesInfo.unidades = ['Gy (total) / Gy-Eq (tumor, normal), ' ...
                'doseScale aplicado'];
        end
    end
    dosesInfo.exportado = true;
    dosesInfo.componentes = listaDosis;
    dosesInfo.doseScale = escalaDosis;
    dosesInfo.archivos = nombresNii;
elseif generarDosis && isempty(doseData)
    warning('BNCT:exportarParaSlicer:sinDoseData', ...
        ['generarDosis=true pero DATA.doseData no existe: no se exporta ' ...
         'ningun .nii (el bundle sigue siendo valido).']);
end

% =========================================================================
% 8) transform.json
% =========================================================================
rutaJson = fullfile(outDir, 'transform.json');
escribirTransformJson(rutaJson, espejarRAS, haySDF, dosesInfo, poseInfo);
nArch = nArch + 1;
archivos{nArch} = rutaJson;
archivos = archivos(1:nArch);                  % recorta la preallocacion

% =========================================================================
% 9) resumen en consola
% =========================================================================
fprintf('\n=== exportarParaSlicer: %d archivos en %s ===\n', ...
    numel(archivos), outDir);
fprintf('  espejarRAS = %d | dosis: %s\n', espejarRAS, ...
    resumenDosis(dosesInfo));
for i = 1:K
    fprintf(['  pose %d: p = [%g %g %g %g %g %g] (mm, grados) | ' ...
        'dmin = %8.3f mm | colisiona: %s\n'], i, P(i, :), ...
        poseInfo(i).dmin_mm, cadenaSiNo(poseInfo(i).colisiona));
end

% =========================================================================
info = struct();
info.outDir = outDir;
info.espejarRAS = espejarRAS;
info.safetyDistance = safetyDistance;
info.nPoses = K;
info.poses = poseInfo;
info.doses = dosesInfo;
info.archivos = archivos;
end

% =========================================================================
% ============================ funciones locales ==========================
% =========================================================================
function ruta = escribirVTK(ruta, pts, caras, titulo)
%ESCRIBIRVTK Escribe un VTK POLYDATA ASCII, formato identico al de
%   data/cone_beam_con_pared_v4_1.vtk: header '# vtk DataFile Version 4.2',
%   titulo, 'ASCII', 'DATASET POLYDATA', POINTS n double y caras como
%   lineas '3 a b c' (POLYGONS m, m*4). Sin caras -> VERTICES (nube de
%   puntos), que Slicer tambien carga como modelo.
n = size(pts, 1);
fid = fopen(ruta, 'w');
if fid == -1
    error('BNCT:exportarParaSlicer:noSePuedeEscribir', ...
        'No se pudo crear %s.', ruta);
end
cleanupVTK = onCleanup(@() fclose(fid));       % cierra incluso con error
fprintf(fid, '# vtk DataFile Version 4.2\n');
fprintf(fid, '%s\n', titulo);
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET POLYDATA\n');
fprintf(fid, 'POINTS %d double\n', n);
fprintf(fid, '%.6g %.6g %.6g\n', pts.');       % un punto por linea
if isempty(caras)
    fprintf(fid, 'VERTICES %d %d\n', n, 2 * n);
    fprintf(fid, '1 %d\n', (1:n).');
else
    m = size(caras, 1);
    fprintf(fid, 'POLYGONS %d %d\n', m, m * 4);   % 3 indices + 1 contador
    fprintf(fid, '3 %d %d %d\n', caras.');        % una cara por linea
end
clear cleanupVTK                               % cierra el archivo
end

% =========================================================================
function conAffine = escribirNIfTI(V, ruta, M, espaciado)
%ESCRIBIRNIFTI Escribe el volumen de dosis y adjunta el affine RAS como
%   sform (Slicer lo interpreta como voxel 0-based -> RAS).
%   En R2023b niftiwrite exige un struct de niftiinfo: se escribe el
%   volumen, se lee el header y se reescribe con TransformName='Sform' y
%   Transform=affine3d(M.') (affine3d usa post-multiplicacion; niftiwrite
%   guarda las filas de M en srow_x/y/z). Si el affine falla, el .nii queda
%   escrito igual y el aviso remite a transform.json (origen/spacing).
conAffine = false;
if exist('niftiwrite', 'file') ~= 2
    warning('BNCT:exportarParaSlicer:sinIPT', ...
        ['niftiwrite no disponible (Image Processing Toolbox): se omite ' ...
         '%s. Instale IPT o use el origen/spacing de transform.json.'], ruta);
    return
end
try
    niftiwrite(single(V), ruta);               % single: universal en Slicer
catch err
    warning('BNCT:exportarParaSlicer:niftiEscritura', ...
        'No se pudo escribir %s (%s).', ruta, err.message);
    return
end
try
    hdr = niftiinfo(ruta);
    hdr.PixelDimensions = abs(espaciado(:).'); % pixdim en mm
    hdr.SpaceUnits = 'Millimeter';
    hdr.TransformName = 'Sform';               % setXForm exige mayuscula inicial
    hdr.Transform = affine3d(M.');             % convencion post-multiplicacion
    niftiwrite(single(V), ruta, hdr);
    conAffine = true;
catch err
    warning('BNCT:exportarParaSlicer:affineFallido', ...
        ['No se pudo escribir el affine sform de %s (%s). El .nii queda ' ...
         'sin posicion: use origen_RAS/spacing_RAS de transform.json o ' ...
         'cargue la dosis con un voxel-to-RAS manual.'], ruta, err.message);
end
end

% =========================================================================
function [D, tipoDosis] = construirDosisCompuesta(nombre, doseData, pesosRBE, escala)
%CONSTRUIRDOSISCOMPUESTA Arma la dosis compuesta en la grilla PHITS.
%   'total'  : suma simple de las 4 componentes fisicas (B+H+N+gamma), en
%              Gy. Es la composicion de calculateUTCP con pesos unitarios
%              (la misma dosis fisica que usa el costo via
%              objectiveFunction). x doseScale (el loader NO lo aplica).
%   'tumor'  : Deq = [D_B D_H D_N D_gamma] * rbeTumor (Gy-Eq), la misma
%              composicion lineal que calculateUTCP hace por voxel.
%   'normal' : idem con rbeNormal.
%   ORDEN DE EJES: los arrays son [nx ny nz] ndgrid (dim1<->x). NIfTI-1
%   guarda dim1 como indice mas rapido y el affine mapea dim1->x: el orden
%   ya coincide con NIfTI. NO permutar.
switch lower(nombre)
    case 'total'
        D = doseData.dose_B + doseData.dose_H + doseData.dose_N + ...
            doseData.dose_gamma;
        tipoDosis = 'Gy';
    case 'tumor'
        w = pesosRBE.tumor(:).';
        D = w(1) * doseData.dose_B + w(2) * doseData.dose_H + ...
            w(3) * doseData.dose_N + w(4) * doseData.dose_gamma;
        tipoDosis = 'Gy-Eq';
    case 'normal'
        w = pesosRBE.normal(:).';
        D = w(1) * doseData.dose_B + w(2) * doseData.dose_H + ...
            w(3) * doseData.dose_N + w(4) * doseData.dose_gamma;
        tipoDosis = 'Gy-Eq';
    otherwise
        error('BNCT:exportarParaSlicer:componenteInvalida', ...
            'Componente de dosis desconocida: %s', nombre);
end
D = D * escala;                                % Gy/source -> Gy (doseScale)
end

% =========================================================================
function [listaDosis, pesosRBE] = planificarDosis(componentes, doseData, ...
    radioBio, generarDosis)
%PLANIFICARDOSIS Resuelve que componentes de dosis se exportan.
%   Regla del proyecto: si se pide 'tumor'/'normal' y radioBio no trae el
%   peso RBE valido (1x4), esa componente se descarta y se asegura 'total'
%   (suma simple x doseScale) con aviso.
listaDosis = {};
pesosRBE = struct('tumor', [], 'normal', []);
if ~generarDosis
    return
end
if isempty(doseData)
    return                                     % el aviso se emite en el flujo
end

pidieronPesos = false;
for c = 1:numel(componentes)
    switch componentes{c}
        case 'total'
            listaDosis{end+1} = 'total';       %#ok<AGROW>
        case 'tumor'
            pidieronPesos = true;
            w = obtenerPesosRBE(radioBio, 'rbeTumor');
            if ~isempty(w)
                pesosRBE.tumor = w;
                listaDosis{end+1} = 'tumor';   %#ok<AGROW>
            end
        case 'normal'
            pidieronPesos = true;
            w = obtenerPesosRBE(radioBio, 'rbeNormal');
            if ~isempty(w)
                pesosRBE.normal = w;
                listaDosis{end+1} = 'normal';  %#ok<AGROW>
            end
        otherwise
            error('BNCT:exportarParaSlicer:componenteInvalida', ...
                ['doseComponentes solo acepta ''total'', ''tumor'' y ' ...
                 '''normal''; recibio ''%s''.'], componentes{c});
    end
end

if pidieronPesos && isempty(pesosRBE.tumor) && isempty(pesosRBE.normal)
    warning('BNCT:exportarParaSlicer:sinPesosRBE', ...
        ['radioBio no trae rbeTumor/rbeNormal validos (1x4): se exporta ' ...
         'solo ''total'' (suma simple de las 4 componentes x doseScale).']);
    if ~any(strcmp(listaDosis, 'total'))
        listaDosis{end+1} = 'total';
    end
end
end

% =========================================================================
function dmin = distanciaExacta(ptsHaz, wallStruct, bodyStruct)
%DISTANCIAEXACTA Minimo de la distancia punto-malla exacta (sin signo, mm)
%   contra wallStruct y bodyStruct (misma estrategia de collisionCheck sin
%   SDF y de validatePHITS).
d = inf(size(ptsHaz, 1), 1);
if mallaValida(wallStruct)
    dw = pointToMeshDistance(ptsHaz, wallStruct.points, wallStruct.faces);
    d = min(d, dw);
end
if mallaValida(bodyStruct)
    db = pointToMeshDistance(ptsHaz, bodyStruct.points, bodyStruct.faces);
    d = min(d, db);
end
dmin = min(d);
end

% =========================================================================
function escribirPoseTxt(ruta, indice, K, p, T, dmin, colisiona, metodo, ...
    safetyDistance)
%ESCRIBIRPOSETXT Pose legible: p, dmin y la matriz T 4x4 (formato %8.3f).
fid = fopen(ruta, 'w');
if fid == -1
    error('BNCT:exportarParaSlicer:noSePuedeEscribir', ...
        'No se pudo crear %s.', ruta);
end
cleanupTxt = onCleanup(@() fclose(fid));
fprintf(fid, 'Pose %d de %d\n', indice, K);
fprintf(fid, ['p = [%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f]   ', ...
    '%% [x y z] mm, [alfa beta gamma] grados\n'], p);
fprintf(fid, 'dmin = %8.3f mm (metodo: %s; seguridad = %g mm)\n', ...
    dmin, metodo, safetyDistance);
fprintf(fid, 'colisiona: %s\n', cadenaSiNo(colisiona));
fprintf(fid, ['T (paciente -> haz, SIN espejo RAS; ultima columna es la ' ...
    'traslacion en mm):\n']);
fprintf(fid, '[%8.3f %8.3f %8.3f %8.3f]\n', T.');
clear cleanupTxt                               % cierra el archivo
end

% =========================================================================
function escribirTransformJson(ruta, espejarRAS, haySDF, dosesInfo, poseInfo)
%ESCRIBIRTRANSFORMJSON transform.json con convencion RAS, grilla de dosis
%   y poses (T SIN espejar, dmin y colision). Todo el texto es ASCII para
%   que fwrite produzca JSON UTF-8 valido sin depender de la codepage.
json = struct();
json.generado = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
json.generador = 'BNCT_PositionOptimizer / exportarParaSlicer';

nota = ['Los VTK exportados (puerto y pacientes) niegan x e y respecto de ' ...
    'los datos de entrada: M_RAS = diag([-1,-1,1,1]) (LPS->RAS). La dosis ' ...
    'NIfTI no se remuestrea: lleva el mismo espejo en el affine sform que ' ...
    'Slicer interpreta como voxel(0-based) -> RAS. T_beam_haz (aqui y en ' ...
    'pose_*.txt) queda SIN espejar: mapea paciente -> haz. ' ...
    'Si tu escena original de Slicer cargaba los .vtk de entrada sin ' ...
    'espejar, regenera con espejarRAS=false.'];
if ~espejarRAS
    nota = ['espejarRAS=false: los archivos quedan en el MARCO DEL HAZ ' ...
        '(como los .vtk de entrada de data/), NO en RAS. ' nota];
end
if ~haySDF
    nota = [nota ' dmin_mm de las poses es distancia EXACTA sin signo ' ...
        '(no hay SDF precomputado): no distingue penetracion (negativa).'];
end
if isfield(dosesInfo, 'conAffine') && ~dosesInfo.conAffine && ...
        dosesInfo.exportado
    nota = [nota ' El .nii de dosis quedo sin affine (ver aviso): use ' ...
        'origen_RAS/spacing_RAS para posicionarlo.'];
end

json.convencion = struct( ...
    'marco_salida', 'RAS', ...
    'espejarRAS', espejarRAS, ...
    'M_RAS_LPSaRAS', 'diag([-1,-1,1,1])', ...
    'nota', nota);

dosesJ = struct('exportado', false, 'origen_RAS', [], 'spacing_RAS', [], ...
    'dims', [], 'unidades', '', 'doseScale', 1, 'componentes', {{}}, ...
    'archivos', {{}});
if dosesInfo.exportado
    dosesJ.exportado = true;
    dosesJ.origen_RAS = dosesInfo.origen_RAS;
    dosesJ.spacing_RAS = dosesInfo.spacing_RAS;
    dosesJ.dims = dosesInfo.dims;
    dosesJ.unidades = dosesInfo.unidades;
    dosesJ.doseScale = dosesInfo.doseScale;
    dosesJ.componentes = dosesInfo.componentes;
    dosesJ.archivos = dosesInfo.archivos;
end
json.doses = dosesJ;

% cell de structs para que jsonencode produzca SIEMPRE un array (con K==1
% un struct array 1x1 se serializaria como objeto unico).
poseCell = cell(1, numel(poseInfo));
for i = 1:numel(poseInfo)
    poseCell{i} = struct('indice', poseInfo(i).indice, ...
        'p', poseInfo(i).p, ...
        'T_beam_haz', poseInfo(i).T, ...
        'dmin_mm', poseInfo(i).dmin_mm, ...
        'colisiona', poseInfo(i).colisiona, ...
        'dminMetodo', poseInfo(i).dminMetodo);
end
json.poses = poseCell;

cuerpo = char(jsonencode(json));
fid = fopen(ruta, 'w');
if fid == -1
    error('BNCT:exportarParaSlicer:noSePuedeEscribir', ...
        'No se pudo crear %s.', ruta);
end
cleanupJson = onCleanup(@() fclose(fid));
fwrite(fid, cuerpo);                           % ASCII -> UTF-8 exacto
clear cleanupJson                              % cierra el archivo
end

% =========================================================================
function P = normalizarPoses(poses)
%NORMALIZARPOSES Acepta Kx6, vector 1x6 o 6x1, cell de 1x6 y struct array
%   con campo .p (mismo criterio que validatePHITS). Devuelve Kx6 double.
if isnumeric(poses)
    P = double(poses);
    if numel(P) == 6
        P = reshape(P, 1, 6);                  % una sola pose (fila o columna)
    end
    if size(P, 2) ~= 6
        error('BNCT:exportarParaSlicer:formatoPoses', ...
            'POSES debe ser Kx6 [x y z alfa beta gamma]; recibio %dx%d.', ...
            size(P, 1), size(P, 2));
    end
elseif iscell(poses)
    P = zeros(numel(poses), 6);
    for i = 1:numel(poses)
        v = poses{i}(:).';
        if numel(v) ~= 6
            error('BNCT:exportarParaSlicer:formatoPoses', ...
                'Cada pose en la cell debe tener 6 elementos; la %d tiene %d.', ...
                i, numel(v));
        end
        P(i, :) = v;
    end
elseif isstruct(poses) && all(isfield(poses, 'p'))
    P = zeros(numel(poses), 6);
    for i = 1:numel(poses)
        P(i, :) = poses(i).p(:).';
    end
else
    error('BNCT:exportarParaSlicer:formatoPoses', ...
        ['POSES debe ser Kx6, un vector 1x6/6x1, una cell de 1x6 o un ' ...
         'struct array con campo .p.']);
end
if any(~isfinite(P), 'all')
    error('BNCT:exportarParaSlicer:poseNoFinita', ...
        'POSES contiene NaN o Inf.');
end
end

% =========================================================================
function componentes = normalizarComponentes(valor)
%NORMALIZARCOMPONENTES char/string/cellstr -> cellstr lowercase sin duplicados.
if ischar(valor)
    valor = {valor};
elseif isstring(valor)
    valor = cellstr(valor);
elseif ~iscell(valor)
    error('BNCT:exportarParaSlicer:componenteInvalida', ...
        'doseComponentes debe ser char, string o cellstr.');
end
componentes = {};
for c = 1:numel(valor)
    nombre = lower(char(valor{c}));
    if ~ismember(nombre, {'total', 'tumor', 'normal'})
        error('BNCT:exportarParaSlicer:componenteInvalida', ...
            ['doseComponentes solo acepta ''total'', ''tumor'' y ' ...
             '''normal''; recibio ''%s''.'], nombre);
    end
    if ~ismember(nombre, componentes)
        componentes{end+1} = nombre;           %#ok<AGROW>
    end
end
end

% =========================================================================
function w = obtenerPesosRBE(params, campo)
%OBTENERPESOSRBE Peso [B H N gamma] de radioBio; [] si falta o es invalido.
%   (calculateUTCP usa [1 1 1 1] como reserva; aqui la reserva es 'total',
%   segun la especificacion del exportador.)
w = [];
if isfield(params, campo) && isnumeric(params.(campo)) && ...
        numel(params.(campo)) == 4
    w = double(params.(campo)(:)).';
end
end

% =========================================================================
function espaciado = espaciadoDeGrilla(gridX, gridY, gridZ)
%ESPACIADODeGRILLA Espaciado (mm) de cada eje; promedio si no es uniforme
%   (con aviso). Rechaza grillas degradadas (<2 puntos o no ascendentes).
grids = {gridX, gridY, gridZ};
nombres = {'gridX', 'gridY', 'gridZ'};
espaciado = zeros(1, 3);
for a = 1:3
    v = grids{a}(:);
    if numel(v) < 2
        error('BNCT:exportarParaSlicer:grillaDegenerada', ...
            '%s necesita al menos 2 puntos (tiene %d).', nombres{a}, numel(v));
    end
    dv = diff(v);
    if any(dv <= 0)
        warning('BNCT:exportarParaSlicer:grillaNoAscendente', ...
            '%s no es ascendente: el affine NIfTI puede quedar inconsistente.', ...
            nombres{a});
    end
    if max(abs(dv - mean(dv))) > 1e-6 * max(1, abs(mean(dv)))
        warning('BNCT:exportarParaSlicer:grillaNoUniforme', ...
            ['%s no es uniforme: se usa el espaciado PROMEDIO (%.4g mm) ' ...
             'en el affine NIfTI.'], nombres{a}, mean(dv));
    end
    espaciado(a) = mean(dv);
end
end

% =========================================================================
function validarDoseData(doseData)
%VALIDARDOSADATA Grillas y arrays de dosis consistentes (ndgrid, igual que
%   sampleDoseMaps/loadPHITSDose).
campos = {'gridX', 'gridY', 'gridZ', 'dose_B', 'dose_H', 'dose_N', 'dose_gamma'};
for c = 1:numel(campos)
    if ~isfield(doseData, campos{c}) || isempty(doseData.(campos{c}))
        error('BNCT:exportarParaSlicer:doseFaltante', ...
            'doseData no tiene el campo requerido "%s".', campos{c});
    end
end
nx = numel(doseData.gridX);
ny = numel(doseData.gridY);
nz = numel(doseData.gridZ);
comps = {'dose_B', 'dose_H', 'dose_N', 'dose_gamma'};
for c = 1:4
    if ~isequal(size(doseData.(comps{c})), [nx ny nz])
        error('BNCT:exportarParaSlicer:doseInconsistente', ...
            '%s es [%s]; la grilla exige [%d %d %d] (ndgrid).', comps{c}, ...
            num2str(size(doseData.(comps{c}))), nx, ny, nz);
    end
end
end

% =========================================================================
function validarCaras(caras, nPts, nombreCampo)
%VALIDARCARAS Mx3 con indices 1-based enteros dentro de [1..nPts].
if size(caras, 2) ~= 3
    error('BNCT:exportarParaSlicer:carasMalFormadas', ...
        '%s.faces debe ser Mx3; recibio %dx%d.', nombreCampo, ...
        size(caras, 1), size(caras, 2));
end
if any(caras(:) < 1) || any(caras(:) > nPts) || ...
        any(caras(:) ~= round(caras(:)))
    error('BNCT:exportarParaSlicer:carasInvalidas', ...
        '%s.faces referencia indices fuera de [1..%d] o no enteros.', ...
        nombreCampo, nPts);
end
end

% =========================================================================
function ok = mallaValida(S)
%MALLAVALIDA Struct con points y faces no vacios (criterio de collisionCheck).
ok = isstruct(S) && isfield(S, 'points') && isfield(S, 'faces') && ...
    ~isempty(S.faces) && ~isempty(S.points);
end

% =========================================================================
function ok = tieneSDF(S)
%TIENESDF Mismo criterio que collisionCheck: sdf + grillas no vacias.
ok = isfield(S, 'sdf') && isfield(S, 'gridX') && isfield(S, 'gridY') && ...
    isfield(S, 'gridZ') && ~isempty(S.sdf);
end

% =========================================================================
function v = valorConDefault(S, campo, def)
%VALORCONDEFAULT Campo de un struct con default (misma idea que
%   objectiveFunction/validatePHITS).
if isfield(S, campo) && ~isempty(S.(campo))
    v = S.(campo);
else
    v = def;
end
end

% =========================================================================
function s = cadenaSiNo(b)
if b
    s = 'si';
else
    s = 'no';
end
end

% =========================================================================
function s = resumenDosis(dosesInfo)
%RESUMENDOSIS Texto corto del estado de la dosis para el resumen final.
if dosesInfo.exportado
    s = strjoin(dosesInfo.componentes, ' + ');
else
    s = '(no exportada)';
end
end
