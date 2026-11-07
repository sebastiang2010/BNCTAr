function [cost, info] = objectiveFunction(p, data)
%OBJECTIVEFUNCTION Funcion de costo del optimizador de posicion (a MINIMIZAR).
%   [cost, info] = OBJECTIVEFUNCTION(P, DATA)
%
%   P    : [x y z alpha beta gamma] (mm, grados); acepta fila o columna.
%   DATA : struct con
%     .patientStruct  : anatomia (ver transformPatient)
%     .wallStruct     : reactor, idealmente con SDF precomputado adjunto
%     .bodyStruct     : reactor adicional opcional (puede faltar)
%     .doseData       : 4 mapas PHITS (marco del haz), ver loadPHITSDose
%     .radioBio       : parametros de calculateUTCP
%     .safetyDistance : mm (default 1)
%     .options        : opcional: .log (container.Map) y/o
%                       .onEvaluation (function handle @(p,cost,info))
%
%   COSTO (escala creciente de gravedad, finito siempre para los solvers):
%     dmin < safetyDistance        -> 1e6 + 1e3*(safety - dmin)   motivo 'colision'
%     mascara tumoral/OAR vacia    -> 1e5                          'tumor_vacio'/'oar_vacio'
%     limite NTCP violado          -> 1e5 + 1e3*exceso             'ntcp_mucosa'/'ntcp_piel'
%     caso aceptable               -> -UTCP                        'ok'
%
%   info : struct con .p .T .dmin .collision .tcp .ntcp .utcp .reason
%          .limitsViolated .detail (salida completa de calculateUTCP).
%
%   REGISTRO DE EVALUACIONES: los structs de MATLAB son por valor, por lo
%   que esta funcion NO puede acumular en data. Dos mecanismos:
%     1) data.options.log como container.Map (clase handle): cada
%        evaluacion agrega la clave 'eval_000001', 'eval_000002', ...
%     2) data.options.onEvaluation: callback @(p, cost, info).
%   optimizePosition ademas envuelve esta funcion en una closure interna
%   que acumula el historial completo (ver documentacion ahi).

if nargin < 2
    error('BNCT:objectiveFunction:faltanArgs', 'Se requieren (p, data).');
end

p = p(:).';
if numel(p) ~= 6
    error('BNCT:objectiveFunction:pInvalido', ...
        'p debe tener 6 elementos [x y z alpha beta gamma]; recibio %d.', numel(p));
end

safetyDistance = 1;   % mm
if isfield(data, 'safetyDistance') && ~isempty(data.safetyDistance)
    safetyDistance = double(data.safetyDistance);
end
siHay = @(S, campo, def) valorConDefault(S, campo, def);

info = struct();
info.p = p;
info.T = rigidTransform(p);
info.reason = 'sin_evaluar';
info.limitsViolated = false;

% ---------------- 1) colision rapida contra el reactor (SDF) ----------------
wallStruct = siHay(data, 'wallStruct', struct('points', zeros(0,3), 'faces', zeros(0,3)));
bodyStruct = siHay(data, 'bodyStruct', struct('points', zeros(0,3), 'faces', zeros(0,3)));
[collision, dmin] = collisionCheck(data.patientStruct, wallStruct, bodyStruct, ...
    info.T, safetyDistance);
info.collision = collision;
info.dmin = dmin;

if collision
    cost = 1e6 + 1e3 * max(0, safetyDistance - dmin);   % penalizacion proporcional
    info.reason = 'colision';
    info.utcp = -Inf;
    info.tcp = NaN;
    info.ntcp = struct('mucosa', NaN, 'piel', NaN);
    info.detail = struct('ok', false, 'reason', 'colision');
    registrar(data, p, cost, info);
    return
end

% ---------------- 2) muestreo de dosis en la anatomia transformada ----------------
[patientT, sampledDose] = transformPatient(data.patientStruct, info.T, data.doseData);

% ---------------- 3) metricas radiobiologicas ----------------
[tcp, ntcp, utcp, detalle] = calculateUTCP(sampledDose, patientT, data.radioBio);
info.tcp = tcp;
info.ntcp = ntcp;
info.utcp = utcp;
info.detail = detalle;

% ---------------- 4) mascaras vacias -> rechazo ----------------
if ~detalle.ok
    cost = 1e5;
    info.reason = normalizarMotivo(detalle.reason);
    registrar(data, p, cost, info);
    return
end

% ---------------- 5) limites clinicos de NTCP ----------------
N_max_m = siHay(data.radioBio, 'N_max_mucosa', 1);
N_max_p = siHay(data.radioBio, 'N_max_piel', 1);
excesoM = max(0, ntcp.mucosa - N_max_m);
excesoP = max(0, ntcp.piel  - N_max_p);

if excesoM > 0 || excesoP > 0
    info.limitsViolated = true;
    cost = 1e5 + 1e3 * max(excesoM, excesoP);
    if excesoM >= excesoP
        info.reason = 'ntcp_mucosa';
    else
        info.reason = 'ntcp_piel';
    end
    registrar(data, p, cost, info);
    return
end

% ---------------- 6) caso aceptable ----------------
info.reason = 'ok';
cost = -utcp;
registrar(data, p, cost, info);
end

% ========================================================================
function registrar(data, p, cost, info)
% Registro por container.Map (handle: persiste fuera de la funcion) y/o
% callback opcional.
if isfield(data, 'options') && isstruct(data.options)
    if isfield(data.options, 'log') && isa(data.options.log, 'container.Map')
        clave = sprintf('eval_%06d', data.options.log.Count + 1);
        data.options.log(clave) = struct('p', p, 'cost', cost, ...
            'reason', info.reason, 'dmin', info.dmin);
    end
    if isfield(data.options, 'onEvaluation') && ...
            isa(data.options.onEvaluation, 'function_handle')
        data.options.onEvaluation(p, cost, info);
    end
end
end

function motivo = normalizarMotivo(reason)
% Reduce motivos compuestos ('oar_vacio:mucosa') al vocabulario de conteo.
if ~ischar(reason)
    motivo = 'otros';
    return
end
pos = strfind(reason, ':');
if ~isempty(pos)
    reason = reason(1:pos(1)-1);
end
switch reason
    case {'tumor_vacio', 'oar_vacio', 'oar_faltante'}
        if startsWith(reason, 'oar')
            motivo = 'oar_vacio';
        else
            motivo = reason;
        end
    otherwise
        motivo = 'otros';
end
end

function v = valorConDefault(S, campo, def)
if isfield(S, campo) && ~isempty(S.(campo))
    v = S.(campo);
else
    v = def;
end
end
