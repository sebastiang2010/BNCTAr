function data = makeDemoData()
%MAKEDEMODATA Genera el set de datos sinteticos completo para el demo.
%   data = MAKEDEMODATA() arma todo lo que objectiveFunction/optimizePosition
%   necesitan, SIN datos reales y sin display:
%
%     - Paciente elipsoidal (semiejes 30x22x25 mm) con tumor esferico
%       excentrico (r=9 mm en (8,9,3)) y 2 OARs como capas 3D:
%         mucosa : capa interna rho in (0.80, 0.92]
%         piel   : capa externa rho in (0.92, 1.00]
%     - Reactor sintetico NO convexo: tubo de pared gruesa (x in [40,90] mm,
%       radio exterior 60 mm, abertura circular de radio 25 mm a lo largo
%       del eje x = eje del haz) como struct points/faces cerrado.
%     - Dosis sintetica en el MARCO DEL HAZ: 4 componentes con atenuacion
%       exponencial axial desde la cara del puerto (x=40, el haz viaja hacia
%       -x) y perfil lateral gaussiano, con constantes distintas por
%       componente. Aproximacion didactica de un mesh tally PHITS.
%     - radioBio: TCP H&N Perotti 2022 y NTCP mucosa Gonzalez 2017 tomados
%       de model_list.txt del proyecto real (../model_list.txt); NTCP piel
%       es un placeholder didactico; pesos RBE/CBE son placeholders
%       (PUNTO DE INTEGRACION, ver README).
%     - safetyDistance = 1 mm y pose inicial de setup pInicial desalineada.
%
%   SALIDA: struct con .patientStruct .wallStruct .bodyStruct .doseData
%           .radioBio .safetyDistance .pInicial .options

res = 2;   % mm (voxel del paciente y de la dosis sintetica)

% ============================ PACIENTE ============================
sx = 30; sy = 22; sz = 25;                        % semiejes (mm)
gx = -sx:res:sx;
gy = -sy:res:sy;
gz = -sz:res:sz;
[X, Y, Z] = ndgrid(gx, gy, gz);
rho = (X/sx).^2 + (Y/sy).^2 + (Z/sz).^2;          % elipsoide normalizado

cuerpo = double(rho <= 1);
tumorC = [8, 9, 3];                               % excentrico: obliga a mover
rTumor = 9;                                       % mm
tumor  = double((X-tumorC(1)).^2 + (Y-tumorC(2)).^2 + (Z-tumorC(3)).^2 <= rTumor^2) .* cuerpo;
mucosa = double(rho > 0.80 & rho <= 0.92) .* cuerpo;
piel   = double(rho > 0.92 & rho <= 1.00);

patientStruct = struct();
patientStruct.mask      = cuerpo;
patientStruct.tumorMask = tumor;
patientStruct.oars      = struct();
patientStruct.oars.mucosa = struct('mask', mucosa, 'nombre', 'mucosa oral');
patientStruct.oars.piel   = struct('mask', piel,   'nombre', 'piel (entrada del haz)');
patientStruct.voxelSize   = [res res res];
patientStruct.origin      = [gx(1), gy(1), gz(1)];
patientStruct.gridX = gx;  patientStruct.gridY = gy;  patientStruct.gridZ = gz;

% superficie del cuerpo (mm, marco paciente) — misma convencion que loadPatient
% (isosurface de 2 salidas devuelve [CARAS, VERTICES] en ese orden)
[faces, verts] = isosurface(cuerpo, 0.5);
patientStruct.surface = struct('points', (double(verts)-1) .* [res res res] ...
    + patientStruct.origin, 'faces', double(faces));

% ============================ REACTOR ============================
% Tubo con abertura circular (NO convexo: la pared interior del orificio es
% una concavidad que alphaShape cerraria). Malla cerrada 2-manifold.
[ptsReactor, fcsReactor] = construirTubo(40, 90, 25, 60, 36);   % x0,x1,rInt,rExt,nAz
wallStruct = struct('points', ptsReactor, 'faces', fcsReactor, ...
    'fuente', 'sintetico: tubo con abertura circular (eje del haz = x)');
bodyStruct = struct('points', zeros(0,3), 'faces', zeros(0,3), ...
    'fuente', '(no provisto)');

% ===================== DOSIS SINTETICA (marco del haz) =====================
% Perfil: D_k(x,y,z) = A_k * exp(-max(0,40-x)/lambda_k) * exp(-(y^2+z^2)/(2*sigma_k^2))
% El haz viaja hacia -x desde la cara del puerto (x=40). Unidades: Gy.
xg = -100:res:100;
yg = -85:res:85;
zg = -85:res:85;
[Xd, Yd, Zd] = ndgrid(xg, yg, zg);
profundidad = max(0, 40 - Xd);   % mm de "tejido equivalente" desde el puerto

comp = @(A, lambda, sigma) A * exp(-profundidad/lambda) .* ...
    exp(-(Yd.^2 + Zd.^2) / (2*sigma^2));
doseB      = comp(3.9, 45, 18);   % boron: componente dominante
doseH      = comp(0.8, 55, 25);   % neutrones rapidos
doseN      = comp(1.2, 40, 20);   % neutrones termicos
doseGamma  = comp(1.5, 60, 22);   % fotones

doseData = loadPHITSDose(struct('gridX', xg, 'gridY', yg, 'gridZ', zg, ...
    'dose_B', doseB, 'dose_H', doseH, 'dose_N', doseN, 'dose_gamma', doseGamma));

% ===================== RADIOBIOLOGIA (demo) =====================
% TCP H&N (Perotti 2022, model_list.txt TCP MODEL 1) y NTCP mucosa
% (Gonzalez 2017, model_list.txt NTCP MODEL 1, a=0 -> Dmax).
% NTCP piel: placeholder didactico (EUD, a=1). Pesos rbe/CBE: placeholder.
radioBio = struct();
radioBio.rbeTumor  = [6.0 3.0 2.0 1.0];   % [B H N gamma] (CBE tumoral alto)
radioBio.rbeNormal = [3.0 3.0 2.0 1.0];   % [B H N gamma] (tejido sano)
radioBio.tcp = struct('Ar', 0.0218, 'Br', 0.00218, 'c1', 2.435, 'c2', 0.120, ...
    'tof', 34.6320, 'tos', 1.2121e+03, ...
    'pf_g', 0.53, 'ps_g', 0.47, 'pf_bnct', 0.2, 'ps_bnct', 0.8);
radioBio.ntcp_mucosa = struct('TD50', 39.8, 'm', 0.17, 'ab', 10, 'a', 0);
radioBio.ntcp_piel   = struct('TD50', 55.0, 'm', 0.50, 'ab', 10, 'a', 1);
radioBio.w_m = 0.4;            % exponente mucosa en UTCP
radioBio.w_p = 0.2;            % exponente piel en UTCP
radioBio.N_max_mucosa = 0.10;  % limite clinico demo
radioBio.N_max_piel   = 0.20;

% ============================ PAQUETE ============================
options = struct();
options.log    = containers.Map;   % registro de evaluaciones (handle)
options.bundleDir = fullfile(tempdir, 'BNCT_validacion_demo');
options.robustDx = 2;              % mm  (robustez validatePHITS)
options.robustDtheta = 2;          % grados

data = struct();
data.patientStruct  = patientStruct;
data.wallStruct     = wallStruct;
data.bodyStruct     = bodyStruct;
data.doseData       = doseData;
data.radioBio       = radioBio;
data.safetyDistance = 1;                                   % mm
data.pInicial       = [10, -8, 6, 8, -6, 4];               % setup desalineado
data.options        = options;
end

% ========================================================================
function [pts, fcs] = construirTubo(x0, x1, rInt, rExt, nAz)
% Malla cerrada de un tubo a lo largo de x: anillos exterior/interior en las
% dos caras (x0 y x1) + cilindro exterior + orificio + 2 coronas anulares.
theta = (0:nAz-1)' * 2*pi / nAz;
cE = [cos(theta), sin(theta)];            % nAz x 2 (y,z)
% vertices: 1..n   exterior@x0 | n+1..2n exterior@x1 | 2n+1..3n interior@x0 | 3n+1..4n interior@x1
vE0 = [repmat(x0, nAz, 1), rExt*cE];
vE1 = [repmat(x1, nAz, 1), rExt*cE];
vI0 = [repmat(x0, nAz, 1), rInt*cE];
vI1 = [repmat(x1, nAz, 1), rInt*cE];
pts = [vE0; vE1; vI0; vI1];

sig = [2:nAz, 1];                          % vecino circular
fcs = zeros(0, 3);
for i = 1:nAz
    j = sig(i);
    % cilindro exterior
    fcs = [fcs; i j nAz+j;  i nAz+j nAz+i];        %#ok<AGROW>
    % orificio (pared interior)
    fcs = [fcs; 2*nAz+i 2*nAz+j 3*nAz+j;  2*nAz+i 3*nAz+j 3*nAz+i]; %#ok<AGROW>
    % corona anular cara x0 (exterior@x0 <-> interior@x0)
    fcs = [fcs; 2*nAz+i 2*nAz+j j;  2*nAz+i j i];  %#ok<AGROW>
    % corona anular cara x1 (exterior@x1 <-> interior@x1)
    fcs = [fcs; nAz+i nAz+j 3*nAz+j;  nAz+i 3*nAz+j 3*nAz+i]; %#ok<AGROW>
end
fcs = double(fcs);
end
