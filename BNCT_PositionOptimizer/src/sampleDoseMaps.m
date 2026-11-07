function dose = sampleDoseMaps(doseData, points)
%SAMPLEDOSEMAPS Muestrea las 4 componentes de dosis en puntos del marco del haz.
%   dose = SAMPLEDOSEMAPS(DOSEDATA, POINTS)
%
%   DOSEDATA : struct (ver loadPHITSDose.m) con:
%     .gridX/.gridY/.gridZ : vectores 1xn en mm, formato NGRID
%                           (gridX varia con la dimension 1, etc.)
%     .dose_B, .dose_H, .dose_N, .dose_gamma : arrays 3D [nx ny nz] en Gy,
%                           formato ndgrid (dim1=x, dim2=y, dim3=z).
%   POINTS   : Nx3 coordenadas YA EN EL MARCO DEL HAZ (mm).
%
%   SALIDA: struct con
%     .B .H .N .gamma : Nx1 dosis interpolada (Gy) en cada punto
%     .nFueraDominio  : cantidad de puntos fuera del dominio de la grilla
%
%   FUERA DE DOMINIO: los puntos fuera de la grilla reciben 0 Gy (interpn
%   con extrapval=0). Es una aproximacion conservadora para el objetivo
%   (menos dosis => menor UTCP) y esta documentada en el README. Cualquier
%   NaN presente en los mapas se propaga como NaN; el limpiado NaN->0 lo
%   hace transformPatient (contando los casos).
%
%   NOTA DE FORMATO: interpn exige formato ndgrid. Usar meshgrid aqui
%   invierte los ejes y produce un muestreo incorrecto (bug corregido).

campos = {'gridX','gridY','gridZ','dose_B','dose_H','dose_N','dose_gamma'};
for i = 1:numel(campos)
    if ~isfield(doseData, campos{i})
        error('BNCT:sampleDoseMaps:faltaCampo', ...
            'doseData no tiene el campo requerido "%s".', campos{i});
    end
end

gX = doseData.gridX(:).';
gY = doseData.gridY(:).';
gZ = doseData.gridZ(:).';
nx = numel(gX); ny = numel(gY); nz = numel(gZ);
esperado = [nx, ny, nz];

componentes = {'dose_B','dose_H','dose_N','dose_gamma'};
for i = 1:numel(componentes)
    if ~isequal(size(doseData.(componentes{i})), esperado)
        error('BNCT:sampleDoseMaps:tamanoinconsistente', ...
            ['"%s" tiene tamano [%s] pero la grilla exige [%d %d %d] ' ...
             '(formato ndgrid).'], componentes{i}, ...
            num2str(size(doseData.(componentes{i}))), nx, ny, nz);
    end
end

q = double(points);
if size(q, 2) ~= 3
    error('BNCT:sampleDoseMaps:puntosMalFormados', ...
        'points debe ser Nx3 (mm); recibido %dx%d.', size(q,1), size(q,2));
end

dose = struct();
dose.B      = interpn(gX, gY, gZ, doseData.dose_B,      q(:,1), q(:,2), q(:,3), 'linear', 0);
dose.H      = interpn(gX, gY, gZ, doseData.dose_H,      q(:,1), q(:,2), q(:,3), 'linear', 0);
dose.N      = interpn(gX, gY, gZ, doseData.dose_N,      q(:,1), q(:,2), q(:,3), 'linear', 0);
dose.gamma  = interpn(gX, gY, gZ, doseData.dose_gamma,  q(:,1), q(:,2), q(:,3), 'linear', 0);

fuera = q(:,1) < gX(1)  | q(:,1) > gX(end) | ...
        q(:,2) < gY(1)  | q(:,2) > gY(end) | ...
        q(:,3) < gZ(1)  | q(:,3) > gZ(end);
dose.nFueraDominio = sum(fuera);
end
