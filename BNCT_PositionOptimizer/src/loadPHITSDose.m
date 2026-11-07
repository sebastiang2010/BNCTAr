function doseData = loadPHITSDose(source)
%LOADPHITSDOSE Carga las 4 matrices de dosis PHITS en el marco del haz.
%   doseData = LOADPHITSDOSE(SOURCE) acepta:
%     (a) ruta a un .mat que contiene un struct (o campos sueltos) con las
%         grillas y componentes (formato "doseData", ver abajo);
%     (b) ruta a un DIRECTORIO con un .mat por componente (nombres que
%         contengan boro/fast/termico/gamma y, opcionalmente, grid);
%     (c) un struct ya armado en memoria (se normaliza y valida).
%
%   FORMATO doseData (CANONICO del proyecto, ver README):
%     .gridX/.gridY/.gridZ : 1xn, mm, ejes NGRID (gridX <-> dimension 1)
%     .dose_B     : [nx ny nz] Gy, dosis por boron   (PHITS: boro)
%     .dose_H     : [nx ny nz] Gy, neutrones rapidos (PHITS: fast)
%     .dose_N     : [nx ny nz] Gy, neutrones termicos (PHITS: thn)
%     .dose_gamma : [nx ny nz] Gy, fotones           (PHITS: Dg / photon)
%     .origin     : [1x3] mm (derivado de la grilla)
%     .spacing    : [1x3] mm (si la grilla es uniforme; [] si no)
%     .orientation: registro de convencion de ejes
%
%   CONVENCIONES DE MARCO: las matrices quedan FIJAS en el MARCO DEL HAZ
%   (estrategia Nivel A). El paciente se mapea al haz con T(p) de
%   rigidTransform; la dosis se muestrea en D_k(T*r).
%
%   INTEGRACION con los parsers del usuario (en el directorio padre):
%     - f_lectura_phits_neutron.m devuelve D1c con .boro/.fast/.thn y
%       .xc/.yc/.zc  ->  mapeo automatico a dose_B/dose_H/dose_N/gridX...
%     - f_lectura_phits_gamma.m devuelve D1c con .Dg  ->  dose_gamma.
%     - lecturaDosis.m (MCNP) devuelve ejeX/ejeY/ejeZ + Matrix 4D con
%       .type: se reconocen los ejes; las componentes deben extraerse a
%       campos con nombres arriba antes de llamar esta funcion.
%     Si SOURCE contiene archivos .out crudos de PHITS se emite un error
%     claro pidiendo convertirlos primero con los parsers del usuario.
%
%   Los NaN de los mapas (regiones de aire sin tally) se reemplazan por
%   0 Gy y se cuenta cuantos fueron (campo .nNaN).

if isstruct(source) && ~isstring(source) && ~ischar(source)
    doseData = normalizar(source);
    return
end
if ~(ischar(source) || isstring(source))
    error('BNCT:loadPHITSDose:sourceInvalido', ...
        'SOURCE debe ser ruta (.mat o directorio) o struct ya armado.');
end
source = char(source);

if exist(source, 'dir') == 7
    doseData = cargarDirectorio(source);
elseif exist(source, 'file') == 2
    doseData = cargarMat(source);
else
    error('BNCT:loadPHITSDose:noExiste', 'Ruta inexistente: %s', source);
end
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
function doseData = cargarDirectorio(dirPath)
archivos = dir(fullfile(dirPath, '*.mat'));
if isempty(archivos)
    outCrudos = dir(fullfile(dirPath, '*.out'));
    if ~isempty(outCrudos)
        error('BNCT:loadPHITSDose:outCrudo', ...
            ['Se encontraron archivos .out de PHITS en "%s". Convierta ' ...
             'primero cada uno con los parsers del proyecto real ' ...
             '(f_lectura_phits_neutron.m para boro/fast/termico y ' ...
             'f_lectura_phits_gamma.m para fotones) y guarde .mat, o ' ...
             'arme el struct doseData a mano.'], dirPath);
    end
    error('BNCT:loadPHITSDose:directorioSinMat', ...
        'No hay archivos .mat de dosis en "%s".', dirPath);
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
