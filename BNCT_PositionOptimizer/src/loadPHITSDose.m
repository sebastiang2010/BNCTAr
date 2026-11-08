function doseData = loadPHITSDose(source, options)
%LOADPHITSDOSE Carga las 4 matrices de dosis PHITS en el marco del haz.
%   doseData = LOADPHITSDOSE(SOURCE)
%   doseData = LOADPHITSDOSE(SOURCE, OPTIONS)
%
%   SOURCE acepta:
%     (0) struct con campos 'neutronFile' y/o 'photonFile' (rutas a .out
%         crudos de PHITS)  ->  MODO REAL: parsea los archivos con los
%         helpers vendidos del usuario:
%           - neutronFile -> leerPHITSNeutron  (adaptado no-interactivo de
%             f_lectura_phits_neutron.m): mset1 -> dose_B (Boron),
%             mset2 -> dose_H (FastNeutron), mset3 -> dose_N (ThermalNeutron).
%           - photonFile  -> leerPHITSGamma    (adaptado de
%             f_lectura_phits_gamma.m): mset1 -> dose_gamma.
%         neutronFile es obligatorio en modo real; si falta photonFile se
%         completa dose_gamma con ceros y se avisa.
%     (a) ruta a un .out crudo de PHITS (neutrones o fotones) -> modo real
%         de un solo archivo (detecta 'part' del header); el resto de las
%         componentes se completa con ceros y se avisa.
%     (b) ruta a un .mat que contiene un struct (o campos sueltos) con las
%         grillas y componentes (formato "doseData", ver abajo);
%     (c) ruta a un DIRECTORIO: si contiene .mat por componente los usa
%         (nombres que contengan boro/fast/termico/gamma y, opcionalmente,
%         grid); si no hay .mat pero hay .out cuyo nombre contenga
%         'neutron' y/o ('photon' o 'gamma'), los carga en modo real.
%     (d) un struct ya armado en memoria (se normaliza y valida).
%
%   FORMATO doseData (CANONICO del proyecto, ver README):
%     .gridX/.gridY/.gridZ : 1xn, mm, ejes NGRID (gridX <-> dimension 1)
%     .dose_B     : [nx ny nz] dosis por boron   (PHITS: mset 1 Boron)
%     .dose_H     : [nx ny nz] neutrones rapidos (PHITS: mset 2 FastNeutron)
%     .dose_N     : [nx ny nz] neutrones termicos(PHITS: mset 3 ThermalNeutron)
%     .dose_gamma : [nx ny nz] fotones           (PHITS: archivo photon)
%     .origin     : [1x3] mm (derivado de la grilla)
%     .spacing    : [1x3] mm (si la grilla es uniforme; [] si no)
%     .orientation: registro de convencion de ejes
%     .doseScale  : factor de escalado DECIDIDO POR EL USUARIO (default 1).
%                   El loader NO lo aplica a las matrices: los .out reales
%                   traen Gy/source (por particula fuente) y la conversion
%                   a dosis absoluta depende del numero de fuentes del
%                   recalculo PHITS. Ver README seccion 6.
%     .units      : 'Gy' (sintetico) o 'Gy/source' (real)
%     .frame      : 'haz' -> las matrices quedan FIJAS en el marco del haz
%                   (marco del tally PHITS con el cilindro; coincide con
%                   el marco de destino de rigidTransform: con p=0 de
%                   rigidTransform las coordenadas de paciente se leen
%                   directamente en este marco).
%     .sourceFiles: cell de archivos de origen (modo real)
%     .meta       : struct(s) devueltos por los helpers (part, dims,
%                   tiempos de parseo, etc.)
%
%   OPTIONS (opcional):
%     .doseScale : factor a registrar en doseData.doseScale (default 1).
%     .verbose   : true/false, imprime resumen (default false).
%
%   CONVENCION DE EJES (ndgrid): las matrices son D(i,j,k) alineadas con
%   gridX(i), gridY(j), gridZ(k), tal que interpn(gridX,gridY,gridZ,D,px,
%   py,pz) interpola correctamente (ver sampleDoseMaps). Los helpers
%   vendidos del usuario ya entregan las paginas en ese orden (x-rapido,
%   y volteada a ascendente segun el header PHITS "y = ny,1,-1"; ver
%   comentarios de leerPHITSNeutron.m / leerPHITSGamma.m).
%
%   INTEGRACION con los parsers del usuario (en el directorio padre):
%     - f_lectura_phits_neutron.m / f_lectura_phits_gamma.m son scripts
%       interactivos (uigetfile); su logica esta vendida en
%       src/leerPHITSNeutron.m y src/leerPHITSGamma.m para uso automatico.
%     - Si SOURCE contiene archivos .out sin nombre reconocible en modo
%       directorio, se emite un error claro.
%
%   Los NaN de los mapas (regiones de aire sin tally) se reemplazan por
%   0 Gy y se cuenta cuantos fueron (campo .nNaN).

if nargin < 2
    options = struct();
end
if ~isstruct(options)
    error('BNCT:loadPHITSDose:opcionesInvalidas', 'OPTIONS debe ser un struct.');
end

% ---------------------------- struct en memoria ----------------------------
if isstruct(source) && ~isstring(source) && ~ischar(source)
    if isfield(source, 'neutronFile') || isfield(source, 'photonFile')
        doseData = cargarPHITSReal(source, options);   % modo real (0)
    else
        doseData = normalizar(source);                 % modo sintetico (d)
    end
    return
end

% ------------------------------ rutas --------------------------------------
if ~(ischar(source) || isstring(source))
    error('BNCT:loadPHITSDose:sourceInvalido', ...
        'SOURCE debe ser ruta (.out/.mat o directorio) o struct ya armado.');
end
source = char(source);

if exist(source, 'file') == 2
    [~,~,ext] = fileparts(source);
    if strcmpi(ext, '.out')
        doseData = cargarPHITSReal(archivarSpec(source), options);  % (a)
    else
        doseData = cargarMat(source);                               % (b)
    end
elseif exist(source, 'dir') == 7
    doseData = cargarDirectorio(source, options);                   % (c)
else
    error('BNCT:loadPHITSDose:noExiste', 'Ruta inexistente: %s', source);
end
end

% ========================================================================
function spec = archivarSpec(rutaOut)
% Envolvemos una ruta .out suelta en el spec del modo real.
nombre = lower(fullfile(rutaOut));
if contains(nombre, 'neutron')
    spec = struct('neutronFile', rutaOut);
elseif contains(nombre, 'photon') || contains(nombre, 'gamma')
    spec = struct('photonFile', rutaOut);
else
    error('BNCT:loadPHITSDose:outSinIdentificar', ...
        ['El archivo "%s" no dice "neutron" ni "photon/gamma" en el nombre. ' ...
         'Use un struct {neutronFile, photonFile} para identificarlo.'], rutaOut);
end
end

% ========================================================================
function doseData = cargarPHITSReal(spec, options)
% MODO REAL: parsea los .out crudos con los helpers vendidos del usuario
% (leerPHITSNeutron / leerPHITSGamma) y arma el struct canonico.
if ~isfield(spec, 'neutronFile')
    error('BNCT:loadPHITSDose:faltaNeutron', ...
        ['El modo real requiere spec.neutronFile (el archivo de neutrones ' ...
         'aporta las 3 componentes B/H/N).']);
end
archivoN = char(spec.neutronFile);
if exist(archivoN, 'file') ~= 2
    error('BNCT:loadPHITSDose:noExiste', 'No existe el archivo de neutrones: %s', archivoN);
end

fprintf('loadPHITSDose (modo real): parseando %s ...\n', archivoN);
[Dn, metaN] = leerPHITSNeutron(archivoN);
fprintf('  neutrones OK: %dx%dx%d, %.1f s (msets 1=Boron, 2=FastNeutron, 3=ThermalNeutron)\n', ...
    size(Dn.boro), metaN.tiempoParseo);

% grilla en mm (los helpers conservan la convencion cm del usuario)
gridX = double(Dn.xc(:)).' * 10;
gridY = double(Dn.yc(:)).' * 10;
gridZ = double(Dn.zc(:)).' * 10;

S = struct('gridX', gridX, 'gridY', gridY, 'gridZ', gridZ, ...
    'dose_B', double(Dn.boro), 'dose_H', double(Dn.fast), 'dose_N', double(Dn.thn));

archivos = {archivoN};
metas = {metaN};

% --- fotones (opcional) ---
if isfield(spec, 'photonFile') && ~isempty(spec.photonFile)
    archivoG = char(spec.photonFile);
    if exist(archivoG, 'file') ~= 2
        error('BNCT:loadPHITSDose:noExiste', 'No existe el archivo de fotones: %s', archivoG);
    end
    fprintf('loadPHITSDose (modo real): parseando %s ...\n', archivoG);
    [Dg, metaG] = leerPHITSGamma(archivoG);
    fprintf('  fotones OK: %dx%dx%d, %.1f s\n', size(Dg.Dg), metaG.tiempoParseo);

    % la grilla del photon debe coincidir con la de neutrones
    if ~isequal(size(Dg.Dg), size(Dn.boro)) || ...
            ~grillasIguales(Dg.xc, Dn.xc) || ~grillasIguales(Dg.yc, Dn.yc) || ...
            ~grillasIguales(Dg.zc, Dn.zc)
        error('BNCT:loadPHITSDose:grillasDiscrepan', ...
            ['La grilla del archivo de fotones (%dx%d en x/y, zdel=%g) no ' ...
             'coincide con la de neutrones (%dx%d en x/y, zdel=%g). Verificar ' ...
             'que ambos tallies usen la misma mesh xyz.'], ...
            numel(Dg.xc), numel(Dg.yc), metaG.xdel, numel(Dn.xc), numel(Dn.yc), metaN.xdel);
    end
    S.dose_gamma = double(Dg.Dg);
    archivos{end+1} = archivoG; %#ok<AGROW>
    metas{end+1} = metaG; %#ok<AGROW>
else
    warning('BNCT:loadPHITSDose:sinPhoton', ...
        ['No se indico photonFile: dose_gamma se completa con 0 ' ...
         '(el calculo de dosis total quedara sin componente gamma).']);
    S.dose_gamma = zeros(size(Dn.boro));
end

% --- campos opcionales del modo real ---
S.doseScale   = obtenerCampo(options, 'doseScale', 1);
S.units       = 'Gy/source';
S.frame       = 'haz';
S.sourceFiles = archivos;
S.meta        = metas;

doseData = normalizar(S);

if obtenerCampo(options, 'verbose', false)
    fprintf(['  grilla (mm): X [%.2f, %.2f] step %.2f | Y [%.2f, %.2f] step %.2f | ' ...
        'Z [%.2f, %.2f] step %.2f\n'], ...
        doseData.gridX(1), doseData.gridX(end), doseData.spacing(1), ...
        doseData.gridY(1), doseData.gridY(end), doseData.spacing(2), ...
        doseData.gridZ(1), doseData.gridZ(end), doseData.spacing(3));
    fprintf('  max por componente [B H N gamma] (Gy/source): %s\n', ...
        mat2str([max(doseData.dose_B(:)), max(doseData.dose_H(:)), ...
                 max(doseData.dose_N(:)), max(doseData.dose_gamma(:))], 4));
end
end

% ========================================================================
function iguales = grillasIguales(a, b)
a = double(a(:)); b = double(b(:));
if numel(a) ~= numel(b)
    iguales = false; return
end
tol = 1e-6 * max(1, max(abs(a)));
iguales = all(abs(a - b) <= tol);
end

% ========================================================================
function doseData = cargarMat(path)
S = load(path);
if isempty(fieldnames(S))
    error('BNCT:loadPHITSDose:matVacio', 'El .mat "%s" esta vacio.', path);
end
% si hay un unico struct, descender en el (caso tipico: D1c o doseData)
campos = fieldnames(S);
if isscalar(S.(campos{1})) && isstruct(S.(campos{1})) && numel(campos) == 1
    S = S.(campos{1});
end
doseData = normalizar(S);
end

% ========================================================================
function doseData = cargarDirectorio(dirPath, options)
% Modo directorio: prioriza .mat por componente; si no hay, intenta modo
% real con los .out de neutrones/fotones presentes.
archivos = dir(fullfile(dirPath, '*.mat'));
if isempty(archivos)
    outNeutron = dir(fullfile(dirPath, '*neutron*.out'));
    outPhoton  = [dir(fullfile(dirPath, '*photon*.out')); dir(fullfile(dirPath, '*gamma*.out'))];
    if ~isempty(outNeutron) || ~isempty(outPhoton)
        % modo real automatico con los .out presentes
        spec = struct();
        if ~isempty(outNeutron)
            spec.neutronFile = fullfile(dirPath, outNeutron(1).name);
        end
        if ~isempty(outPhoton)
            spec.photonFile = fullfile(dirPath, outPhoton(1).name);
        end
        if ~isfield(spec, 'neutronFile')
            warning('BNCT:loadPHITSDose:soloPhoton', ...
                ['El directorio "%s" solo tiene .out de fotones: dose_B/H/N ' ...
                 'quedaran en 0.'], dirPath);
        end
        doseData = cargarPHITSReal(spec, options);
        return
    end
    error('BNCT:loadPHITSDose:directorioSinMat', ...
        'No hay archivos .mat ni .out reconocibles (neutron/photon) en "%s".', dirPath);
end

S = struct();
for i = 1:numel(archivos)
    parcial = cargarMat(fullfile(dirPath, archivos(i).name));
    campos = fieldnames(parcial);
    for c = 1:numel(campos)
        destino = mapear(lower(campos{c}));
        if ~isempty(destino) && ~isfield(S, destino)
            S.(destino) = parcial.(campos{c});
        end
    end
end
doseData = normalizar(S);
end

% ========================================================================
function doseData = normalizar(S)
S = aplanar(S);   % si S tiene un unico campo struct, descender
salida = struct();

requeridos = {'gridX','gridY','gridZ','dose_B','dose_H','dose_N','dose_gamma'};
faltantes = '';
camposOrigen = fieldnames(S);
for r = 1:numel(requeridos)
    if isfield(S, requeridos{r})
        salida.(requeridos{r}) = double(S.(requeridos{r}));
    else
        encontrado = false;
        for c = 1:numel(camposOrigen)
            clave = mapear(lower(camposOrigen{c}));
            if strcmp(clave, requeridos{r})
                salida.(requeridos{r}) = double(S.(camposOrigen{c}));
                encontrado = true;
                break
            end
        end
        if ~encontrado
            faltantes = [faltantes, ' ', requeridos{r}]; %#ok<AGROW>
        end
    end
end
if ~isempty(faltantes)
    error('BNCT:loadPHITSDose:camposFaltantes', ...
        ['Faltan componentes/grillas del formato doseData:%s.\n' ...
         'Campos detectados: %s\n' ...
         'Nombres alternativos reconocidos: xc/ejeX->gridX, boro->dose_B, ' ...
         'fast->dose_H, thn/termico->dose_N, Dg/gamma->dose_gamma.'], ...
        faltantes, strjoin(camposOrigen, ', '));
end

% --- tamanos consistentes con la grilla (ndgrid) ---
nx = numel(salida.gridX); ny = numel(salida.gridY); nz = numel(salida.gridZ);
comps = {'dose_B','dose_H','dose_N','dose_gamma'};
for i = 1:4
    if ~isequal(size(salida.(comps{i})), [nx ny nz])
        error('BNCT:loadPHITSDose:tamanoinconsistente', ...
            '%s es [%s]; la grilla exige [%d %d %d] (ndgrid).', comps{i}, ...
            num2str(size(salida.(comps{i}))), nx, ny, nz);
    end
end

% --- NaN -> 0 (aire sin tally), contados ---
nNaN = 0;
for i = 1:4
    nNaN = nNaN + sum(isnan(salida.(comps{i})), 'all');
    salida.(comps{i})(isnan(salida.(comps{i}))) = 0;
end

salida.origin      = [salida.gridX(1), salida.gridY(1), salida.gridZ(1)];
difs = [diffUniforme(salida.gridX), diffUniforme(salida.gridY), diffUniforme(salida.gridZ)];
salida.spacing     = difs;
salida.orientation = 'ndgrid-XYZ (dim1=x, dim2=y, dim3=z); marco del haz';
salida.nNaN        = nNaN;

% --- campos descriptivos (modo real y defaults para el sintetico) ---
salida.doseScale   = obtenerCampo(S, 'doseScale', 1);
if isfield(S, 'units') && ~isempty(S.units)
    salida.units = char(S.units);
else
    salida.units = 'Gy';
end
if isfield(S, 'frame') && ~isempty(S.frame)
    salida.frame = char(S.frame);
else
    salida.frame = 'haz';
end
if isfield(S, 'sourceFiles')
    salida.sourceFiles = S.sourceFiles;
else
    salida.sourceFiles = {};
end
if isfield(S, 'meta')
    salida.meta = S.meta;
else
    salida.meta = {};
end

doseData = salida;
end

function destino = mapear(nombreCampoLower)
% Devuelve el nombre canonico si el campo es reconocido; '' si no.
tabla = {'gridx','gridX'; 'xc','gridX'; 'ejex','gridX'; 'x','gridX'; ...
         'gridy','gridY'; 'yc','gridY'; 'ejey','gridY'; 'y','gridY'; ...
         'gridz','gridZ'; 'zc','gridZ'; 'ejez','gridZ'; 'z','gridZ'; ...
         'dose_b','dose_B'; 'boro','dose_B'; 'db','dose_B'; ...
         'dose_h','dose_H'; 'fast','dose_H'; 'rapido','dose_H'; ...
         'dose_n','dose_N'; 'thn','dose_N'; 'termico','dose_N'; ...
         'thermal','dose_N'; 'epitermico','dose_N'; ...
         'dose_gamma','dose_gamma'; 'dg','dose_gamma'; 'gamma','dose_gamma'; ...
         'foton','dose_gamma'; 'photon','dose_gamma'};
destino = '';
fila = strcmp(nombreCampoLower, tabla(:,1));
if any(fila)
    destino = tabla{find(fila, 1), 2};
end
end

function s = aplanar(S)
% Si S es un envoltorio con un unico campo struct, descender un nivel.
campos = fieldnames(S);
if numel(campos) == 1 && isstruct(S.(campos{1})) && isscalar(S.(campos{1}))
    s = S.(campos{1});
else
    s = S;
end
end

function d = diffUniforme(v)
dv = diff(v(:));
if isempty(dv)
    d = 0;
    return
end
if max(abs(dv - dv(1))) > 1e-6 * max(1, abs(dv(1)))
    warning('BNCT:loadPHITSDose:grillaNoUniforme', ...
        'La grilla no es uniforme; spacing guarda solo el primer paso.');
end
d = dv(1);
end

function valor = obtenerCampo(S, campo, def)
if isfield(S, campo) && ~isempty(S.(campo))
    valor = S.(campo);
else
    valor = def;
end
end
