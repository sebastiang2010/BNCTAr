function [patientT, sampledDose] = transformPatient(patientStruct, T, doseData)
%TRANSFORMPATIENT Transforma el paciente a T y muestrea la dosis en sus mascaras.
%   [patientT, sampledDose] = TRANSFORMPATIENT(PATIENTSTRUCT, T, DOSEDATA)
%
%   ESTRATEGIA NIVEL A: los 4 mapas PHITS quedan FIJOS en el MARCO DEL HAZ;
%   para cada T = rigidTransform(p) se muestrea D_k^p(r) = D_k^phantom(T*r).
%   T mapea paciente -> haz: r_haz = R*r_paciente + t.
%
%   PATIENTSTRUCT : con .tumorMask, .oars.(nombre).mask, .surface.points,
%                   .voxelSize, .origin y .gridX/.gridY/.gridZ (mm, ndgrid).
%                   Cada ROI puede traer su propia grilla (.gridX/Y/Z y
%                   .voxelSize dentro de la entrada); si no, se usa la del
%                   cuerpo. Las mascaras son 3D en formato ndgrid
%                   (mask(i,j,k), i<->x).
%   T             : 4x4 de rigidTransform (paciente -> haz).
%   DOSEDATA      : struct canonico de loadPHITSDose.
%
%   SALIDA:
%     patientT    : copia de patientStruct con .surface.points transformado
%                   (para colisiones) y .T adjunta.
%     sampledDose : struct con
%       .tumor               : Nt x 4 [D_B D_H D_N D_gamma] Gy (por voxel)
%       .tumorVoxelVolumeCm3 : escalar
%       .oars.(nombre)       : struct con .dose (Nx4), .mean (1x4),
%                              .voxelVolumeCm3, .voxelsFueraDominio
%       .voxelsFueraDominio.tumor : cantidad de voxels fuera de la grilla
%
%   FUERA DE DOMINIO: voxels que caen fuera de la grilla de dosis reciben
%   0 Gy (interpn extrapval=0, NaN->0). Es la aproximacion documentada del
%   Nivel A (el mapa no extiende la dosis simulada fuera de su dominio).

R = T(1:3, 1:3);
t = T(1:3, 4);

% ---------------- superficie transformada (para colisiones) ----------------
patientT = patientStruct;
patientT.surface.points = (R * patientStruct.surface.points.' + t).';
patientT.T = T;

% ---------------- muestreo de dosis en mascara tumoral ----------------
if ~isfield(patientStruct, 'tumorMask')
    error('BNCT:transformPatient:sinTumor', ...
        ['patientStruct no tiene .tumorMask. Asignelo desde el NIfTI de ' ...
         'la ROI tumoral (ver loadPatient.m y README).']);
end

sampledDose = struct();
sampledDose.oars = struct();
sampledDose.voxelsFueraDominio = struct();

[salidaT, dvT, nFueraT] = muestrearRegion(patientStruct.tumorMask, ...
    obtenerGrillas(patientStruct, 'tumor'), R, t, doseData);
sampledDose.tumor = salidaT;
sampledDose.tumorVoxelVolumeCm3 = dvT;
sampledDose.voxelsFueraDominio.tumor = nFueraT;

% ---------------- muestreo por OAR ----------------
if isfield(patientStruct, 'oars') && ~isempty(fieldnames(patientStruct.oars))
    nombres = fieldnames(patientStruct.oars);
    for i = 1:numel(nombres)
        entrada = patientStruct.oars.(nombres{i});
        if isstruct(entrada) && isfield(entrada, 'mask')
            mascara = entrada.mask;
        else
            mascara = entrada;   % tolera oars.mucosa = matriz directa
        end
        [salidaO, dvO, nFueraO] = muestrearRegion(mascara, ...
            obtenerGrillas(patientStruct, nombres{i}), R, t, doseData);
        reg = struct();
        reg.dose = salidaO;
        if isempty(salidaO)
            reg.mean = zeros(1, 4);
        else
            reg.mean = mean(salidaO, 1);
        end
        reg.voxelVolumeCm3 = dvO;
        reg.voxelsFueraDominio = nFueraO;
        sampledDose.oars.(nombres{i}) = reg;
    end
end
end

% ========================================================================
function [D, dvCm3, nFuera] = muestrearRegion(mascara, grids, R, t, doseData)
% Extrae los centros de voxel de la mascara, los transforma al haz y
% muestrea las 4 componentes.
gx = grids.gridX(:);   % forzar columna: indexar un vector fila conserva
gy = grids.gridY(:);   % su orientacion y romperia la concatenacion Nx3
gz = grids.gridZ(:);
voxelSize = grids.voxelSize;

lin = find(double(mascara) > 0);
if isempty(lin)
    D = zeros(0, 4);
    dvCm3 = prod(double(voxelSize)) / 1000;
    nFuera = 0;
    return
end
% find en 3D solo devuelve 2 subindices + valores: convertir con ind2sub
[ii, jj, kk] = ind2sub(size(mascara), lin);
pts = [gx(ii), gy(jj), gz(kk)];          % Nx3 mm, marco paciente
q = (R * pts.' + t).';                   % Nx3 mm, marco haz

comp = sampleDoseMaps(doseData, q);      % usa interpn con ndgrid (correcto)
D = [comp.B, comp.H, comp.N, comp.gamma];
nFuera = comp.nFueraDominio;
D(~isfinite(D)) = 0;                     % NaN -> 0 (documentado)
dvCm3 = prod(double(voxelSize)) / 1000;  % mm^3 -> cm^3
% nFuera > 0 indica voxels fuera del dominio simulado (dosis 0 Gy); el
% paciente deberia permanecer dentro del dominio para que el Nivel A sea
% confiable. No detiene el calculo.
end

% ========================================================================
function grids = obtenerGrillas(patientStruct, nombreRegion)
% Resuelve la grilla de una region: propia si existe, si no la del cuerpo.
g = struct('gridX', patientStruct.gridX, 'gridY', patientStruct.gridY, ...
           'gridZ', patientStruct.gridZ, 'voxelSize', patientStruct.voxelSize);

if strcmp(nombreRegion, 'tumor') && isfield(patientStruct, 'tumorGrid') ...
        && isstruct(patientStruct.tumorGrid)
    g = completar(g, patientStruct.tumorGrid);
elseif ~strcmp(nombreRegion, 'tumor') && isfield(patientStruct, 'oars') ...
        && isfield(patientStruct.oars, nombreRegion) ...
        && isstruct(patientStruct.oars.(nombreRegion))
    g = completar(g, patientStruct.oars.(nombreRegion));
end
grids = g;
end

function g = completar(g, fuente)
if isfield(fuente, 'gridX'); g.gridX = fuente.gridX; end
if isfield(fuente, 'gridY'); g.gridY = fuente.gridY; end
if isfield(fuente, 'gridZ'); g.gridZ = fuente.gridZ; end
if isfield(fuente, 'voxelSize'); g.voxelSize = fuente.voxelSize; end
end
