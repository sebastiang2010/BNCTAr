%%RUN_OPTIMIZER Demo end-to-end del optimizador de posicion de paciente BNCT.
% Corre SIN datos reales y SIN display (headless-safe). Flujo:
%   1) addpath src y datos sinteticos (makeDemoData)
%   2) SDF del reactor (signedDistanceMap3D) para colisiones rapidas
%   3) linea de base con la pose de setup pInicial (desalineada)
%   4) optimizacion (surrogateopt, ~150 evaluaciones, minutos)
%   5) validacion exacta + robustez de los 3 mejores finalistas
%      (validatePHITS: pointToMeshDistance + bundle para Nivel B)
%   6) resumen impreso en consola

% --------------------------------------------------------------------------
rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'src'));
oldDir = cd(rootDir);                       % cwd estable (el bundle va a options.bundleDir)
restaurar = onCleanup(@() cd(oldDir));

fprintf('\n=== OPTIMIZADOR DE POSICION BNCT (demo sintetica) ===\n');

%% 1) datos sinteticos
data = makeDemoData();
fprintf('Datos sinteticos listos: paciente %dx%dx%d voxels, reactor %d triangulos.\n', ...
    size(data.patientStruct.mask), size(data.wallStruct.faces, 1));

%% 2) SDF del reactor (una sola vez; colisiones rapidas en el bucle)
data.wallStruct = signedDistanceMap3D(data.wallStruct, ...
    struct('resolution', 2, 'margin', [100 60 60]));
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
    % reserva: usar la mejor pose del historial aunque haya sido rechazada
    [~, j] = min(history.cost);
    finalistas = history.p(j, :);
end
end
