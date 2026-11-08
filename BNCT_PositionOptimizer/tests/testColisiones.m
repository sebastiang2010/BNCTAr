%TESTCOLISIONES Validacion del detector de colisiones con casos conocidos.
%   Script sin entrada: construye mallas sinteticas con separaciones
%   conocidas y compara SDF y distancia exacta contra valores de
%   referencia. Separa:
%     A) caja CERRADA (fillMethod='parity'): separacion/contacto/penetracion
%     B) parche ABIERTO (fillMethod='shell'): plano, a ambos lados
%     C) cono REAL: punto critico p=0 (cuerpo tocando pared)
%
%   Criterio de paso: |SDF - referencia| <= 1.5 voxel (banda + interpolacion)
%   y |exacta - referencia| <= 0.05 mm para distancias unsigned.

addpath(fileparts(mfilename('fullpath')));   % .../src
fprintf('=== VALIDACION DETECTOR DE COLISIONES ===\n');
paso = 0; fallo = 0;
verifica = @(nombre, obtenido, ref, tol) verificar(nombre, obtenido, ref, tol);

%% A) Caja cerrada 100 mm (parity) ---------------------------------------
caja = cajaCerrada(100);                        % cubo [-50,50]^3
sdfCaja = signedDistanceMap3D(caja, struct('resolution', 2, 'margin', [40 40 40]));

casosCaja = { % {vector [x y z], referencia esperada (mm), nombre}
    {[0 0 55],  5, 'A1 separacion +5 (afuera)'}
    {[0 0 50],  0, 'A2 contacto (cara z=+50)'}
    {[0 0 47], -3, 'A3 penetracion -3 (adentro)'}
    {[0 0 25], -25, 'A4 interior profundo -25'}
};
fprintf('\n[A] Caja cerrada (parity), voxel 2 mm:\n');
for k = 1:numel(casosCaja)
    pto = casosCaja{k}{1}; ref = casosCaja{k}{2}; nombre = casosCaja{k}{3};
    dSdf  = evaluaSDFinline(sdfCaja, pto);
    [dEx, ~] = pointToMeshDistance(pto, caja.points, caja.faces);
    paso = paso + verifica([nombre ' (SDF)'],    dSdf, ref, 3.0);
    paso = paso + verifica([nombre ' (exacta)'], dEx, abs(ref), 0.05);
end

%% B) Parche abierto: plano cuadrado 200x200 en z=0 (shell) ---------------
plano = planoAbierto(200);
sdfPlano = signedDistanceMap3D(plano, struct('resolution', 2, ...
    'margin', [40 40 40], 'fillMethod', 'shell'));

casosPlano = { % el "negativo" en un parche abierto = lado opuesto de la pared
    {[0 0  4],  4, 'B1 lado exterior +4'}
    {[0 0 -2], -2, 'B2 lado interior -2'}
    {[0 0  1],  1, 'B3 contacto cercano +1 (SDF conservador)'}
};
fprintf('\n[B] Parche abierto (shell), voxel 2 mm:\n');
for k = 1:numel(casosPlano)
    pto = casosPlano{k}{1}; ref = casosPlano{k}{2}; nombre = casosPlano{k}{3};
    dSdf  = evaluaSDFinline(sdfPlano, pto);
    [dEx, ~] = pointToMeshDistance(pto, plano.points, plano.faces);
    paso = paso + verifica([nombre ' (exacta)'], dEx, abs(ref), 0.05);
    if k < 3
        paso = paso + verifica([nombre ' (SDF)'], dSdf, ref, 3.0);
    else
        % B3: punto a 1 mm de la pared = colision con seguridad 1 mm ->
        % el SDF del shell debe ser CONSERVADOR (dSdf < 1), no exacto.
        okCons = dSdf < 1;
        fprintf('  [%s] %-42s SDF=%7.2f (<1 requerido: conservador)\n', ...
            ternario(okCons), nombre, dSdf);
        paso = paso + double(okCons);
    end
end

%% C) Cono real: punto critico del BODY en p=0 ---------------------------
raiz = fileparts(fileparts(mfilename('fullpath')));
conePath = fullfile(raiz, 'data', 'cone_beam_con_pared_v4_1.vtk');
bodyPath = fullfile(raiz, 'data', 'BODY.vtk');
if exist(conePath, 'file') && exist(bodyPath, 'file')
    cono = loadReactorSTL(conePath);
    nuc  = clipMeshBox(cono, [-420 20 -200 200 -200 200]);
    sdfN = signedDistanceMap3D(nuc, struct('resolution', 2.5, ...
        'margin', [60 60 60], 'fillMethod', 'shell'));
    body = loadReactorSTL(bodyPath);

    % punto mas profundo del BODY (sin transformar) segun el SDF
    dTodos = evaluaSDFtodos(sdfN, body.points);
    [dMinSdf, j] = min(dTodos);
    ptoCritico = body.points(j, :);
    [dMinEx, ~] = pointToMeshDistance(ptoCritico, nuc.points, nuc.faces);

    fprintf('\n[C] Cono real (shell, 2.5 mm):\n');
    fprintf('    punto critico BODY: (%.1f, %.1f, %.1f) mm\n', ptoCritico);
    fprintf('    SDF = %.2f mm | exacta unsigned = %.2f mm | dif = %.2f mm\n', ...
        dMinSdf, dMinEx, dMinSdf - dMinEx);
    % semantica: el SDF en modo shell da el lado de la pared con sesgo de
    % banda (~1.5 voxel). La exacta unsigned solo afirma "distancia minima".
    if dMinSdf < 0 && dMinEx < 1.5
        fprintf('    [C1] consistente: SDF dice penetracion/roce, exacta ~0\n');
        paso = paso + 1;
    else
        fprintf('    [C1] INCONSISTENTE: revisar (SDF=%.2f, exacta=%.2f)\n', ...
            dMinSdf, dMinEx);
        fallo = fallo + 1;
    end
else
    fprintf('\n[C] OMITIDA: faltan data/cone_beam_con_pared_v4_1.vtk o BODY.vtk\n');
end

%% Resumen ----------------------------------------------------------------
fprintf('\n=== RESUMEN: %d verificaciones OK, %d fallos ===\n', paso, fallo);
if fallo > 0
    error('BNCT:testColisiones:fallos', '%d casos fallaron', fallo);
end

% ========================================================================
function ok = verificar(nombre, obtenido, ref, tol)
ok = abs(obtenido - ref) <= tol;
if ok
    fprintf('  [OK]   %-42s obtenido=%7.2f esperado=%7.2f\n', nombre, obtenido, ref);
else
    fprintf('  [FALLO]%-42s obtenido=%7.2f esperado=%7.2f (tol=%.2f)\n', ...
        nombre, obtenido, ref, tol);
end
end

% ------------------------------------------------------------------------
function s = ternario(ok)
if ok, s = 'OK'; else, s = 'FALLO'; end
end

% ------------------------------------------------------------------------
function S = cajaCerrada(L)
% Cubo [-L/2, L/2]^3 con 12 triangulos.
h = L/2;
v = [-h -h -h; h -h -h; h h -h; -h h -h; -h -h h; h -h h; h h h; -h h h];
f = [1 3 2; 1 4 3;   % z=-h (cara inferior)
     5 6 7; 5 7 8;   % z=+h (superior)
     1 2 6; 1 6 5;   % y=-h
     3 4 8; 3 8 7;   % y=+h
     2 3 7; 2 7 6;   % x=+h
     4 1 5; 4 5 8];  % x=-h
S = struct('points', v, 'faces', f);
end

% ------------------------------------------------------------------------
function S = planoAbierto(L)
% Plano z=0, dos triangulos (parche abierto).
v = [-L -L 0; L -L 0; L L 0; -L L 0];
f = [1 2 3; 1 3 4];
S = struct('points', v, 'faces', f);
end

% ------------------------------------------------------------------------
function d = evaluaSDFinline(S, pto)
% Valor del SDF en UN punto (interpolacion trilineal).
d = interpn(S.gridX, S.gridY, S.gridZ, S.sdf, pto(1), pto(2), pto(3), 'linear');
end

% ------------------------------------------------------------------------
function d = evaluaSDFtodos(S, P)
% SDF evaluado en N puntos (solo dominio; fuera -> +Inf).
d = inf(size(P, 1), 1);
in = P(:,1) >= S.gridX(1) & P(:,1) <= S.gridX(end) & ...
     P(:,2) >= S.gridY(1) & P(:,2) <= S.gridY(end) & ...
     P(:,3) >= S.gridZ(1) & P(:,3) <= S.gridZ(end);
d(in) = interpn(S.gridX, S.gridY, S.gridZ, S.sdf, ...
    P(in,1), P(in,2), P(in,3), 'linear');
end
