function [result, history] = optimizePosition(data, opts)
%OPTIMIZEPOSITION Driver de optimizacion de la pose rigida del paciente BNCT.
%   [result, history] = OPTIMIZEPOSITION(DATA, OPTS)
%
%   OPTS (struct opcional):
%     .method    : 'surrogateopt' (default) | 'ga' | 'patternsearch' |
%                  'bayesopt' | 'grid'
%     .lb        : [1x6] limites inferiores (mm, grados)
%                  default [-30 -30 -30 -15 -15 -15]
%     .ub        : [1x6] limites superiores default [30 30 30 15 15 15]
%     .maxFunEvals : presupuesto de evaluaciones (default 500)
%     .x0        : semilla para patternsearch/grid (default centro de lb/ub)
%     .verbose   : true -> Display 'iter' (default false, headless-safe)
%     .plot      : true -> grafico de convergencia si hay display (opcional)
%
%   RESULT:
%     .pBest .TBest .utcp .tcp .ntcp .dmin .cost .counts (por motivo)
%     .nEvals .method .elapsed .lb .ub
%   HISTORY: struct con .p (Nx6) .cost (Nx1) .reason (Nx1 cell) .dmin (Nx1)
%            (una fila por evaluacion).
%
%   El objetivo es la closure @(p) objectiveFunction(p, data) envuelta en
%   una funcion anidada que acumula history (los structs de MATLAB son por
%   valor; la funcion anidada comparte su workspace con este driver).
%   La funcion objetivo JAMAIS devuelve -Inf/NaN: los rechazos son
%   penalizaciones finitas grandes (ver objectiveFunction).

if nargin < 2 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'method') || isempty(opts.method)
    opts.method = 'surrogateopt';
end
if ~isfield(opts, 'lb') || isempty(opts.lb)
    opts.lb = [-30 -30 -30 -15 -15 -15];
end
if ~isfield(opts, 'ub') || isempty(opts.ub)
    opts.ub = [30 30 30 15 15 15];
end
if ~isfield(opts, 'maxFunEvals') || isempty(opts.maxFunEvals)
    opts.maxFunEvals = 500;
end
if ~isfield(opts, 'verbose') || isempty(opts.verbose)
    opts.verbose = false;
end
if ~isfield(opts, 'x0') || isempty(opts.x0)
    opts.x0 = (opts.lb + opts.ub) / 2;
end
if ~isfield(opts, 'plot') || isempty(opts.plot)
    opts.plot = false;
end

lb = double(opts.lb(:)).';
ub = double(opts.ub(:)).';
if any(ub <= lb)
    error('BNCT:optimizePosition:limitesInvalidos', ...
        'Se requiere lb < ub en las 6 componentes.');
end
validarData(data);

displayLevel = 'off';
if opts.verbose
    displayLevel = 'iter';
end

% ---------------- closure con historial (funciones anidadas) ----------------
history = struct('p', zeros(0,6), 'cost', zeros(0,1), 'reason', {{}}, 'dmin', zeros(0,1));
nEval = 0;
tInicio = tic;

    function [c, info] = wrapped(p)
        nEval = nEval + 1;
        [c, info] = objectiveFunction(p, data);
        history.p(end+1, :) = p(:).';
        history.cost(end+1, 1) = c;
        history.reason{end+1, 1} = info.reason;
        history.dmin(end+1, 1) = info.dmin;
    end

    function [c, constraints] = bayesWrap(x)
        p = deBayes(x);
        [c, info] = wrapped(p);
        colisionViol = max(0, data.safetyDistance - info.dmin);
        ntcpViol = 0;
        if info.limitsViolated
            ntcpViol = max(0, info.ntcp.mucosa - data.radioBio.N_max_mucosa);
            ntcpViol = max(ntcpViol, max(0, info.ntcp.piel - data.radioBio.N_max_piel));
        end
        constraints = [colisionViol, ntcpViol];
    end

    function p = deBayes(x)
        % bayesopt entrega x como struct/tabla de variables con nombre
        if istable(x)
            x = table2struct(x);
        end
        if isstruct(x)
            nombres = {'tx','ty','tz','alpha','beta','gamma'};
            p = zeros(1, 6);
            for i = 1:6
                p(i) = double(x.(nombres{i}));
            end
        else
            p = double(x(:).');
        end
    end

    function pBest = gridSearch()
        % Busqueda por barridos coordenados sobre grilla (determinista):
        % pasada 1 en todo el rango, pasada 2 en una ventana reducida.
        pBest = min(max(double(opts.x0(:).'), lb), ub);
        [cBest, ~] = wrapped(pBest);
        nPuntos = 5;
        for pasada = 1:2
            for d = 1:6
                if pasada == 1
                    centro = (lb(d) + ub(d)) / 2;
                    span = (ub(d) - lb(d));
                else
                    centro = pBest(d);
                    span = (ub(d) - lb(d)) / 2;
                end
                vals = centro + linspace(-span/2, span/2, nPuntos);
                vals = min(max(vals, lb(d)), ub(d));
                for v = vals
                    q = pBest;
                    q(d) = v;
                    [c, ~] = wrapped(q);
                    if c < cBest
                        cBest = c;
                        pBest = q;
                    end
                end
            end
        end
    end

metodo = lower(opts.method);

switch metodo
    case 'surrogateopt'
        verificarToolbox('surrogateopt', 'Global Optimization Toolbox');
        options = optimoptions('surrogateopt', ...
            'MaxFunctionEvaluations', opts.maxFunEvals, ...
            'Display', displayLevel, 'PlotFcn', {});
        [pBest, ~] = surrogateopt(@(x) wrapped(x), lb, ub, options);

    case 'ga'
        verificarToolbox('ga', 'Global Optimization Toolbox');
        popSize = 50;
        maxGen = max(5, ceil(opts.maxFunEvals / popSize));
        options = optimoptions('ga', 'PopulationSize', popSize, ...
            'MaxGenerations', maxGen, 'Display', displayLevel, 'PlotFcn', {});
        [pBest, ~] = ga(@(x) wrapped(x), 6, [], [], [], [], lb, ub, [], options);

    case 'patternsearch'
        verificarToolbox('patternsearch', 'Global Optimization Toolbox');
        x0 = min(max(double(opts.x0(:).'), lb), ub);
        options = optimoptions('patternsearch', ...
            'MaxFunctionEvaluations', opts.maxFunEvals, ...
            'Display', displayLevel, 'PlotFcn', {});
        [pBest, ~] = patternsearch(@(x) wrapped(x), x0, [], [], [], [], ...
            lb, ub, [], options);

    case 'bayesopt'
        verificarToolbox('bayesopt', 'Global Optimization Toolbox / Statistics and ML');
        vars = [optimizableVariable('tx',    [lb(1) ub(1)])
                optimizableVariable('ty',    [lb(2) ub(2)])
                optimizableVariable('tz',    [lb(3) ub(3)])
                optimizableVariable('alpha', [lb(4) ub(4)])
                optimizableVariable('beta',  [lb(5) ub(5)])
                optimizableVariable('gamma', [lb(6) ub(6)])];
        resultados = bayesopt(@bayesWrap, vars, ...
            'MaxObjectiveEvaluations', opts.maxFunEvals, ...
            'NumCoupledConstraints', 2, ...
            'IsObjectiveDeterministic', true, ...
            'Verbose', double(opts.verbose), 'PlotFcn', {});
        % Se elige la mejor evaluacion REAL del historial propio en lugar de
        % bestPoint(): con restricciones acopladas bestPoint puede fallar o
        % devolver tipos que deBayes no resuelve (historial siempre Nx6).
        [~, jMin] = min(history.cost);
        pBest = history.p(jMin, :);
        if isempty(pBest)
            error('BNCT:optimizePosition:bayesSinEvaluaciones', ...
                'bayesopt no registro ninguna evaluacion.');
        end

    case 'grid'
        pBest = gridSearch();

    otherwise
        error('BNCT:optimizePosition:metodoDesconocido', ...
            'opts.method "%s" no reconocido.', opts.method);
end

% ---------------- metricas finales consistentes en pBest ----------------
elapsed = toc(tInicio);
[costFinal, infoFinal] = objectiveFunction(pBest, data);

result = struct();
result.pBest   = pBest(:).';
result.TBest   = rigidTransform(result.pBest);
result.utcp    = infoFinal.utcp;
result.tcp     = infoFinal.tcp;
result.ntcp    = infoFinal.ntcp;
result.dmin    = infoFinal.dmin;
result.cost    = costFinal;
result.reason  = infoFinal.reason;
result.counts  = contarMotivos(history.reason);
result.nEvals  = nEval;
result.method  = metodo;
result.elapsed = elapsed;
result.lb      = lb;
result.ub      = ub;

% grafico opcional de convergencia (solo si hay display disponible)
if opts.plot
    try
        figure('Name', 'Convergencia optimizePosition');
        plot(history.cost, 'b.-');
        grid on;
        xlabel('Evaluacion');
        ylabel('Costo (menor es mejor; -UTCP en casos validos)');
        title(sprintf('Metodo: %s | UTCP final = %.4f', metodo, result.utcp));
    catch errGrafico
        warning('BNCT:optimizePosition:sinGrafico', ...
            'No se pudo graficar (entorno sin display): %s', errGrafico.message);
    end
end
end

% ========================================================================
function validarData(data)
campos = {'patientStruct', 'doseData', 'radioBio'};
for i = 1:numel(campos)
    if ~isfield(data, campos{i})
        error('BNCT:optimizePosition:faltaData', ...
            'data no tiene el campo requerido "%s".', campos{i});
    end
end
if ~isfield(data.patientStruct, 'tumorMask')
    error('BNCT:optimizePosition:sinTumor', ...
        'data.patientStruct no tiene .tumorMask (ver loadPatient/README).');
end
if ~isfield(data.patientStruct, 'surface') || ~isfield(data.patientStruct.surface, 'points')
    error('BNCT:optimizePosition:sinSuperficie', ...
        'data.patientStruct.surface.points es requerido para colisiones.');
end
if ~isfield(data, 'safetyDistance') || isempty(data.safetyDistance)
    data.safetyDistance = 1;   % no muta al llamador; objectiveFunction re-default
end
end

% ========================================================================
function counts = contarMotivos(reasons)
motivos = {'ok', 'colision', 'ntcp_mucosa', 'ntcp_piel', 'tumor_vacio', 'oar_vacio', 'otros'};
counts = struct();
sumaConocidos = 0;
for i = 1:numel(motivos) - 1
    counts.(motivos{i}) = sum(strcmp(reasons, motivos{i}));
    sumaConocidos = sumaConocidos + counts.(motivos{i});
end
counts.otros = max(0, numel(reasons) - sumaConocidos);
end

% ========================================================================
function verificarToolbox(funcion, toolbox)
if exist(funcion, 'file') ~= 2 && exist(funcion, 'builtin') ~= 5
    error('BNCT:optimizePosition:sinToolbox', ...
        'El metodo requiere %s (funcion %s no encontrada).', toolbox, funcion);
end
end
