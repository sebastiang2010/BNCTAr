%%RUN_OPTIMIZER Demo end-to-end del optimizador de posicion de paciente BNCT.
% Corre SIN display (headless-safe). Flujo:
%   1) addpath src y datos (makeDemoData +, opcionalmente, dosis real PHITS)
%   2) SDF del reactor (signedDistanceMap3D) para colisiones rapidas
%   3) linea de base con la pose de setup pInicial (desalineada)
%   4) optimizacion (surrogateopt, ~150 evaluaciones, minutos)
%   5) validacion exacta + robustez de los 3 mejores finalistas
%      (validatePHITS: pointToMeshDistance + bundle para Nivel B)
%   6) resumen impreso en consola
%
% MODO DE DATOS: cambie modoDatos a 'real' para usar la dosimetria real
% de cilindro (data/mesh-neutron.out + data/mesh-photon.out, Gy/source)
% cargada con loadPHITSDose, y la geometria real de colision:
%   - cono del reactor: data/cone_beam_con_pared_v4_1.vtk, filtrado con
%     clipMeshBox para descartar la envolvente lejana (basura del modelo).
%   - paciente: la superficie de data/BODY.vtk se usa como
%     patientStruct.surface.points (se mueve con T(p)).
% Las mascaras anatomicas (tumor/OARs) siguen siendo las sinteticas de
% makeDemoData hasta contar con NIfTI reales: UTCP demo-escala, colision
% real. Ver README.md seccion 6: en modo real hay que decidir doseScale
% (los .out traen Gy/source) y los pesos radioBio del demo son placeholders.

% --------------------------------------------------------------------------
modoDatos = 'real';    % 'demo' (default, dosis sintetica) | 'real' (PHITS .out)

rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'src'));
oldDir = cd(rootDir);                       % cwd estable (el bundle va a options.bundleDir)
restaurar = onCleanup(@() cd(oldDir));

fprintf('\n=== OPTIMIZADOR DE POSICION BNCT (modo de datos: %s) ===\n', modoDatos);

%% 1) datos
data = makeDemoData();
fprintf('Datos base listos: paciente %dx%dx%d voxels, reactor %d triangulos.\n', ...
    size(data.patientStruct.mask), size(data.wallStruct.faces, 1));

if strcmpi(modoDatos, 'real')
    % ---- dosimetria real de cilindro (Gy/source) ----
    archivoNeutron = fullfile(rootDir, 'data', 'mesh-neutron.out');
    archivoPhoton  = fullfile(rootDir, 'data', 'mesh-photon.out');
    doseReal = loadPHITSDose(struct('neutronFile', archivoNeutron, ...
                                    'photonFile', archivoPhoton), ...
                             struct('verbose', true));
    data.doseData = doseReal;

    % ---- geometria real de colision (mismo marco PHITS) ----
    % Cono del reactor: se filtra la envolvente lejana (X hasta +-2600 mm,
    % basura del modelo); queda el nucleo del haz X in [-400, 0] mm,
    % radio +-150 mm, que termina donde empieza la grilla PHITS (x = 0).
    coneReal = loadReactorSTL(fullfile(rootDir, 'data', ...
        'cone_beam_con_pared_v4_1.vtk'));
    cajaCono = [-420 20 -200 200 -200 200];   % núcleo + margen de holgura
    data.wallStruct = clipMeshBox(coneReal, cajaCono);

    % Superficie del paciente (BODY.vtk, ~500x250x310 mm): se mueve con T(p)
    % y colisiona contra el cono estatico.
    bodyReal = loadReactorSTL(fullfile(rootDir, 'data', 'BODY.vtk'));
    % .faces no interviene en colisionCheck (usa .points) pero permite que
    % exportarParaSlicer genere la malla del paciente en vez de una nube de
    % puntos (BODY.vtk trae 81139 triangulos).
    data.patientStruct.surface = struct('points', bodyReal.points, 'faces', bodyReal.faces);
    data.bodyStruct = struct('points', zeros(0,3), 'faces', zeros(0,3), ...
        'fuente', '(no provisto en modo real: BODY.vtk es la superficie del paciente)');

    fprintf(['Geometria real: cono %d tri (de %d; envolvente filtrada), ' ...
        'paciente BODY %d pts / %d tri.\n'], ...
        size(data.wallStruct.faces, 1), size(coneReal.faces, 1), ...
        size(bodyReal.points, 1), size(bodyReal.faces, 1));

    % aviso de cobertura: la grilla real debe contener al paciente (tras
    % T(p)); fuera de dominio la dosis interpolada es 0 (interpn extrapval=0)
    gR = [doseReal.gridX(1) doseReal.gridX(end); ...
          doseReal.gridY(1) doseReal.gridY(end); ...
          doseReal.gridZ(1) doseReal.gridZ(end)];
    vox = centrosDeMascara(data.patientStruct.tumorMask, ...
        data.patientStruct.gridX, data.patientStruct.gridY, data.patientStruct.gridZ);
    dentro = vox >= gR(:,1)' & vox <= gR(:,2)';
    fraccionTumor = mean(all(dentro, 2));
    fprintf(['\nAVISO (modo real): la grilla PHITS cubre X [%.1f, %.1f] mm, ' ...
        'Y/Z [%.1f, %.1f] mm. Solo el %.1f%% de los voxels del tumor demo ' ...
        'cae dentro del dominio; el resto recibe 0 Gy (fuera de dominio).\n' ...
        'Las mascaras anatomicas siguen siendo sinteticas (UTCP demo-escala; ' ...
        'colision con geometria real).\n' ...
        'Las dosis crudas son Gy/source (doseScale=%.3g): ajuste doseScale y ' ...
        'radioBio antes de usar metricas clinicas.\n'], ...
        gR(1,1), gR(1,2), gR(2,1), gR(2,2), 100*fraccionTumor, doseReal.doseScale);
    sdfRes = 2.5;  sdfMargin = [60 60 60];    % cono real: caja chica
    % El nucleo del cono es un parche ABIERTO (la boca en x=0 no tiene tapa):
    % fill 'parity' clasificaria el tunel del haz como "interior solido" y
    % marcaria colision fantasma para el paciente dentro del tunel. 'shell'
    % usa rasterizado FINO (solo voxels a <= res/2 del plano de algun
    % triangulo, + dilatacion 1 voxel): el sdf vale ~0 en la pared y crece
    % linealmente hacia afuera; hacia adentro de la banda puede dar ~-6 mm
    % (sesgo conservador del shell, ver signedDistanceMap3D).
    sdfFill = 'shell';
elseif strcmpi(modoDatos, 'demo')
    fprintf('Datos sinteticos listos (dosis sintetica en Gy).\n');
    sdfRes = 2;    sdfMargin = [100 60 60];   % reactor demo grande
    sdfFill = 'parity';                       % tubo demo cerrado (2-manifold)
else
    error('modoDatos invalido: use ''demo'' o ''real''.');
end

%% 2) SDF del reactor (una sola vez; colisiones rapidas en el bucle)
data.wallStruct = signedDistanceMap3D(data.wallStruct, ...
    struct('resolution', sdfRes, 'margin', sdfMargin, 'fillMethod', sdfFill));
fprintf('SDF del reactor: grilla %dx%dx%d @ %g mm (fill=%s).\n', ...
    size(data.wallStruct.sdf), data.wallStruct.resolution, data.wallStruct.fillMethod);

%% 3) linea de base: pose de setup inicial (desalineada)
[~, info0] = objectiveFunction(data.pInicial, data);
fprintf('\n--- Pose de setup inicial p0 = [%g %g %g %g %g %g] ---\n', data.pInicial);
imprimirMetricas(info0);

%% 4) optimizacion
opts = struct();
opts.method = 'surrogateopt';
opts.maxFunEvals = 150;
opts.verbose = false;
[result, history] = optimizePosition(data, opts);

fprintf('\n--- Optimizacion (%s, %d evaluaciones, %.1f s) ---\n', ...
    result.method, result.nEvals, result.elapsed);
fprintf('Mejor pose p* = [%.2f %.2f %.2f %.2f %.2f %.2f]  (mm, grados)\n', result.pBest);
fprintf('Matriz T* (paciente -> haz):\n');
disp(result.TBest);
imprimirMetricas(result);
fprintf('Evaluaciones por motivo: ok=%d, colision=%d, ntcp_mucosa=%d, ntcp_piel=%d, ...\n', ...
    result.counts.ok, result.counts.colision, result.counts.ntcp_mucosa, ...
    result.counts.ntcp_piel);

%% 5) finalistas (top-3 distintos entre si, solo evaluaciones aceptadas)
finalistas = seleccionarFinalistas(history, 3);
fprintf('\n--- Validacion de %d finalistas (colision exacta + robustez) ---\n', ...
    size(finalistas, 1));
report = validatePHITS(finalistas, data);

%% 6) resumen de validacion
for k = 1:size(finalistas, 1)
    fprintf(['Finalista %d: UTCP=%.4f | TCP=%.4f | NTCP mucosa=%.4f piel=%.4f | ' ...
        'dmin exacta=%.2f mm | UTCP perturbado [%.4f, %.4f] | poses limite violadas=%d/12\n'], ...
        report(k).id, report(k).utcp, report(k).tcp, report(k).ntcp.mucosa, ...
        report(k).ntcp.piel, report(k).dminExacta, ...
        report(k).utcpPertMin, report(k).utcpPertMax, report(k).nPertViolan);
end
fprintf('\nBundle para recalculo PHITS (Nivel B): %s\n', report(1).bundleDir);
fprintf('=== FIN DEL DEMO ===\n\n');

%% ---------------- funciones locales ----------------
function vox = centrosDeMascara(mascara, gx, gy, gz)
% Centros de voxel (Nx3, mm) de los voxels con mascara > 0 (marco paciente).
gx = gx(:); gy = gy(:); gz = gz(:);
lin = find(double(mascara) > 0);
if isempty(lin)
    vox = zeros(0, 3);
    return
end
[ii, jj, kk] = ind2sub(size(mascara), lin);
vox = [gx(ii), gy(jj), gz(kk)];
end

function imprimirMetricas(info)
% info: salida de objectiveFunction o struct result de optimizePosition
% (ambos comparten los campos utcp/tcp/ntcp/dmin/reason).
utcp = info.utcp; tcp = info.tcp; ntcp = info.ntcp; dmin = info.dmin;
motivo = info.reason;
if isinf(utcp) && utcp < 0
    estado = 'RECHAZADA';
else
    estado = 'aceptada';
end
fprintf(['  UTCP=%.4f | TCP=%.4f | NTCP mucosa=%.4f | NTCP piel=%.4f | ' ...
    'dmin=%.2f mm | [%s: %s]\n'], utcp, tcp, ntcp.mucosa, ntcp.piel, dmin, estado, motivo);
end

function finalistas = seleccionarFinalistas(history, nDeseados)
% Top-n de poses aceptadas mutuamente separadas (>3 mm o >1.5 grados).
idx = find(strcmp(history.reason, 'ok'));
[~, orden] = sort(history.cost(idx));
idx = idx(orden);
finalistas = zeros(0, 6);
for i = 1:numel(idx)
    p = history.p(idx(i), :);
    if isempty(finalistas)
        separado = true;
    else
        dTr = vecnorm(finalistas(:, 1:3) - p(1:3), 2, 2);
        dRot = vecnorm(finalistas(:, 4:6) - p(4:6), 2, 2);
        separado = all(dTr > 3 | dRot > 1.5);   % separada de TODAS las ya elegidas
    end
    if separado
        finalistas(end+1, :) = p;   %#ok<AGROW>
    end
    if size(finalistas, 1) == nDeseados
        break
    end
end
if isempty(finalistas)
    error(['BNCT:run_optimizer:sinPosicionesFactibles', ...
        ' No se encontraron posiciones factibles (sin colision y con NTCP ' ...
        'dentro de los limites) en %d evaluaciones. Verifique el rango de ' ...
        'busqueda (lb/ub), la geometria (marco de mallas) y la dosis antes ' ...
        'de maximizar UTCP: primero debe existir una region factible.'], nEval);
end
end
