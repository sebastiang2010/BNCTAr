function [tcp, ntcp, utcp, info] = calculateUTCP(dose, patientStruct, params)
%CALCULATEUTCP Calcula TCP, NTCP (mucosa y piel) y UTCP = TCP*(1-NTCPm)^wm*(1-NTCPp)^wp.
%   [tcp, ntcp, utcp, info] = CALCULATEUTCP(DOSE, PATIENTSTRUCT, PARAMS)
%
%   DOSE          : struct muestreado por transformPatient:
%       .tumor                 : Nt x 4 [D_B D_H D_N D_gamma] Gy por voxel
%       .tumorVoxelVolumeCm3   : volumen de voxel del tumor (cm^3)
%       .oars.(nombre)         : struct con .dose (Nx4 Gy), .voxelVolumeCm3
%       Si falta .tumorVoxelVolumeCm3 se deduce de patientStruct.voxelSize.
%   PATIENTSTRUCT : para volumen de voxel de reserva y registro.
%   PARAMS        : struct de radiobiologia ("radioBio"), campos:
%       .rbeTumor  [1x4] pesos Gy->Gy-Eq del tumor  ([B H N gamma])
%       .rbeNormal [1x4] pesos Gy->Gy-Eq de los OARs
%       .tcp         : struct (modelo Perotti 2022, port de TCPc.m):
%                      Ar, Br, c1, c2, tof, tos, pf_g, ps_g, pf_bnct, ps_bnct
%       .tcpValue    : (alternativa) TCP fijo; tiene precedencia si existe
%       .ntcp_mucosa : struct LKB (TD50, m, ab, a)  [port de ntcp.m]
%                      o valor numerico directo de NTCP
%       .ntcp_piel   : idem
%       .w_m, .w_p   : exponentes de penalizacion (adimensional)
%       .N_max_mucosa, .N_max_piel : limites clinicos de NTCP (0-1)
%
%   SALIDAS:
%     tcp  : probabilidad de control tumoral (0-1)
%     ntcp : struct con campos .mucosa y .piel (0-1)
%     utcp : TCP*(1-NTCPm)^w_m * (1-NTCPp)^w_p; = -Inf si un NTCP supera su
%            N_max o si una mascara requerida esta vacia.
%     info : struct con .ok (calculo posible), .reason (motivo de rechazo),
%            .limitsViolated, .deqTumorMean (Gy-Eq), .tumorVoxels.
%
%   MODELOS (ports fieles de los archivos reales del usuario, en ../):
%     - TCP  : TCPc.m  (G F Perotti Bernardini et al 2022 Radiat. Res. 198),
%              con parametros de model_list.txt (TCP H&N).
%     - NTCP : ntcp.m  (S J Gonzalez et al 2017 PMB 62 7938 con a=0 -> Dmax;
%              A M Dattoli Viegas et al 2024 Phys. Med. 128 104840 con a~=0
%              -> EUD).
%
%   PUNTO DE INTEGRACION PROVENZANO (ver README):
%     - La composicion Deq = [D_B D_H D_N D_gamma] * rbe es una suma
%       ponderada SIMPLE. El modelo exacto del usuario (ntcp2.m) combina las
%       componentes con factores de reparacion G y terminos A,B (alpha/beta
%       por componente). Reemplace la seccion marcada
%       "COMPOSICION DEQ (INTEGRACION)" por esa conversion, o llame
%       directamente a ntcp2.m / TCPc.m con los vóxeles muestreados.
%
%   UNIDADES: dosis en Gy (componentes) y Gy-Eq (tras el peso RBE);
%   volumenes en cm^3; probabilidades adimensionales.

% ------------------------- defaults seguros -------------------------
info = struct('ok', false, 'reason', 'sin_evaluar', 'limitsViolated', false, ...
    'deqTumorMean', NaN, 'tumorVoxels', 0);
tcp  = 0;
ntcp = struct('mucosa', NaN, 'piel', NaN);
utcp = -Inf;

rbeT = obtenerRBE(params, 'rbeTumor');
rbeN = obtenerRBE(params, 'rbeNormal');

% ==================== TCP ====================
Dtum = dose.tumor;                         % Nt x 4 Gy
if isempty(Dtum)
    info.reason = 'tumor_vacio';
    return
end
info.tumorVoxels = size(Dtum, 1);

dvT = volumenVoxel(dose, 'tumorVoxelVolumeCm3', patientStruct);   % cm^3

% ----- COMPOSICION DEQ (INTEGRACION): suma ponderada simple.
%       Reemplazable por ntcp2.m (terminos A,B con reparacion) ---
DeqT = Dtum * rbeT(:);                     % Nt x 1 Gy-Eq
info.deqTumorMean = mean(DeqT);

if esCampoNumerico(params, 'tcpValue')
    tcp = double(params.tcpValue);
elseif esCampoStruct(params, 'tcp')
    tcp = tcpPerotti(DeqT, dvT, params.tcp);
else
    % reserva didactica: logistico simple sobre dosis media
    tcp = 1 / (1 + exp(-(0.35 * mean(DeqT) - 4)));
    warning('BNCT:calculateUTCP:tcpGenerico', ...
        'params no trae .tcp ni .tcpValue: se usa un logistico generico.');
end

% ==================== NTCP ====================
% --- mucosa ---
[ntcp.mucosa, motivoM] = evaluaOAR(dose, patientStruct, 'mucosa', rbeN, ...
    params.ntcp_mucosa);
if ~isempty(motivoM)
    info.reason = motivoM;
    return
end

% --- piel ---
[ntcp.piel, motivoP] = evaluaOAR(dose, patientStruct, 'piel', rbeN, ...
    params.ntcp_piel);
if ~isempty(motivoP)
    info.reason = motivoP;
    return
end

% ==================== UTCP ====================
w_m = obtenerEscalar(params, 'w_m', 0);    % exponente mucosa
w_p = obtenerEscalar(params, 'w_p', 0);    % exponente piel
N_max_m = obtenerEscalar(params, 'N_max_mucosa', 1);
N_max_p = obtenerEscalar(params, 'N_max_piel', 1);

utcp = tcp * (1 - ntcp.mucosa)^w_m * (1 - ntcp.piel)^w_p;
info.ok = true;
info.reason = 'ok';

% ----- limites clinicos explicitos: UTCP = -Inf si se violan -----
if ntcp.mucosa > N_max_m || ntcp.piel > N_max_p
    utcp = -Inf;
    info.limitsViolated = true;
    if ntcp.mucosa > N_max_m
        info.reason = 'limite_ntcp:mucosa';
    else
        info.reason = 'limite_ntcp:piel';
    end
end
end

% ========================================================================
function tcp = tcpPerotti(Deq, dv, prm)
% Port fiel de TCPc.m (Perotti Bernardini et al 2022, Radiat. Res. 198).
% Deq: dosis tumoral por voxel en unidades equivalentes a fotones (Gy-Eq);
% dv : volumen de voxel (cm^3).
time = 30;                                  % min (irradiacion de referencia)
xf = prm.tof ./ time;
xs = prm.tos ./ time;
Gf = 2*xf .* (1 - xf .* (1 - exp(-1 ./ xf)));
Gs = 2*xs .* (1 - xs .* (1 - exp(-1 ./ xs)));
% factor de tiempo de referencia de fotones (f_g=1, f_bnct=0, como TCPc.m)
Gr = Gs - (1 * prm.pf_g + 0 * prm.pf_bnct) * (Gs - Gf);

Deq = Deq(~isnan(Deq));
H = exp(-prm.Ar .* Deq - Gr .* prm.Br .* (Deq.^2)).^(1 / prm.c2);
tcp = exp(-prm.c1 * (sum(H) * dv)^prm.c2);
end

% ========================================================================
function [y, motivo] = evaluaOAR(dose, patientStruct, nombre, rbeN, prmNTCP)
% NTCP de un OAR: struct LKB (port de ntcp.m) o valor directo numerico.
motivo = '';
if ~isfield(dose, 'oars') || ~isfield(dose.oars, nombre)
    motivo = ['oar_faltante:' nombre];      % OAR no definido -> fallo seguro
    y = NaN;
    return
end
reg = dose.oars.(nombre);
Deq = reg.dose * rbeN(:);                  % Nx1 Gy-Eq
if isempty(Deq)
    motivo = ['oar_vacio:' nombre];         % mascara sin voxels muestreados
    y = NaN;
    return
end

if isnumeric(prmNTCP)
    y = double(prmNTCP);                    % NTCP provisto directamente
    return
end
if ~isstruct(prmNTCP)
    error('BNCT:calculateUTCP:ntcpParamInvalido', ...
        'params.ntcp_%s debe ser struct (TD50,m,ab,a) o numerico.', nombre);
end

dv = volumenVoxel(reg, 'voxelVolumeCm3', patientStruct);   % cm^3
Deq = Deq(~isnan(Deq));
if all(Deq <= 0)
    y = 0;                                  % sin dosis -> sin complicacion
    return
end

if prmNTCP.a == 0
    % modelo de dosis MAXIMA (mucositis, Gonzalez 2017)
    D = max(Deq);
    EQD2 = D .* (prmNTCP.ab + D) ./ (prmNTCP.ab + 2);
    s = (EQD2 - prmNTCP.TD50) ./ (prmNTCP.m * prmNTCP.TD50);
else
    % modelo EUD (Dattoli 2024)
    vi = dv ./ (dv * numel(Deq));           % fraccion de volumen por voxel
    EUD = (sum(Deq.^prmNTCP.a .* vi))^(1 / prmNTCP.a);
    s = (EUD - prmNTCP.TD50) ./ (prmNTCP.m * prmNTCP.TD50) / 2^0.5;
end
y = 1 ./ 2 * (1 + erf(s));
end

% ========================================================================
function rbe = obtenerRBE(params, campo)
if isfield(params, campo) && numel(params.(campo)) == 4
    rbe = double(params.(campo)(:)).';
else
    rbe = [1 1 1 1];    % sin pesos: Deq = suma de componentes
end
end

function v = obtenerEscalar(params, campo, def)
if isfield(params, campo) && isnumeric(params.(campo)) && ~isempty(params.(campo))
    v = double(params.(campo)(1));
else
    v = def;
end
end

function ok = esCampoNumerico(S, campo)
ok = isfield(S, campo) && isnumeric(S.(campo)) && isscalar(S.(campo));
end

function ok = esCampoStruct(S, campo)
ok = isfield(S, campo) && isstruct(S.(campo));
end

function dv = volumenVoxel(reg, campo, patientStruct)
% Volumen de voxel en cm^3 (mm^3/1000). Reserva: patientStruct.voxelSize.
if isfield(reg, campo) && ~isempty(reg.(campo))
    dv = double(reg.(campo));
elseif isfield(patientStruct, 'voxelSize') && numel(patientStruct.voxelSize) == 3
    dv = prod(double(patientStruct.voxelSize)) / 1000;
else
    dv = 1;   % ultimo recurso: 1 cm^3 (documentado)
end
end
