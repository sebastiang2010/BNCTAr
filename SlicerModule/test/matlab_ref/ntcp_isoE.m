function NTCP = ntcp_isoE(Xiso, Yiso, TD50, m, a, ab)
% ntcp_isoE - NTCP LKB sobre DVH isoefectivo en Gy(IsoE).
%   NTCP = ntcp_isoE(Xiso, Yiso, TD50, m, a, ab) calcula la probabilidad
%   de complicacion con el modelo Lyman-Kutcher-Burman a partir del DVH
%   acumulativo isoefectivo (Xiso en Gy(IsoE), Yiso en %).
%
%   Entradas:
%     Xiso - vector dosis isoefectiva [Gy(IsoE)] (formato dvh_isoE: 0 antepuesto).
%     Yiso - vector volumen acumulativo [%] (mismo largo que Xiso).
%     TD50 - dosis de tolerancia al 50% [Gy(IsoE)] del endpoint.
%     m    - pendiente del modelo LKB.
%     a    - exponente EUD (a=0 -> Dmax; a=1 -> Dmean; a=10 -> criterio cerebral).
%     ab   - cociente alpha/beta fotonico [Gy] (solo rama a=0, conversion EQD2).
%
%   Endpoints registrados (ver md/isoefectiva_B2_B3.md):
%     - Mucosa PROVISIONAL: TD50=39.8, m=0.17, a=0, ab=10
%       (mucositis >= 3, Gonzalez 2017 / Strigari 2012). NO validado para
%       piel, hueso, tronco ni medula: uso provisional y reversible.
%     - Cerebro: TD50=13.9049, m=0.1766, a=10, ab=3
%       (somnolencia SS >= 2, Dattoli 2024).
%
%   Definiciones:
%
%   EUD = (sum_i v_i * D_i^a)^(1/a),  a > 0
%
%   EQD2 = Dmax * (ab + Dmax) / (ab + 2),  a = 0
%
%   t = (Dref - TD50) / (m * TD50),  NTCP = 0.5 * (1 + erf(t / sqrt(2)))
%
%   donde Dref = EUD (a > 0) o EQD2 (a = 0).
%
%   Notas:
%     - La entrada DEBE ser Gy(IsoE) (salida de isodose via recalc_isoE_ROI
%       + dvh_isoE). NO pasar Gy(RBE): el pipeline UTCP v25/v8 usa esta
%       funcion y ya no alimenta el NTCP con dosis RBE.
%     - Desvio documentado respecto de C:\programas\BNCTAr\BNCTAr\ntcp.m
%       (solo lectura): la rama mucosa (a=0) de ntcp.m usa y=0.5*(1+erf(s))
%       sin /sqrt(2), mientras su rama cerebral si incluye /sqrt(2) en s.
%       Aqui se usa erf(t/sqrt(2)) en AMBAS ramas (forma LKB estandar).
%     - NO usar mat/ntcp2.m (incompleto: a_ij sin asignar, alpha_r y tiempo
%       fijos, erf sin /sqrt(2)).
%
%   Ver tambien: dvh_isoE, recalc_isoE_ROI, isoE_params.

if nargin ~= 6
    error('ntcp_isoE: se requieren 6 argumentos (Xiso, Yiso, TD50, m, a, ab).');
end
Xiso = Xiso(:);
Yiso = Yiso(:);
if numel(Xiso) ~= numel(Yiso) || numel(Xiso) < 2
    error('ntcp_isoE: Xiso e Yiso deben tener igual largo (>= 2).');
end
if any(~isfinite(Xiso)) || any(~isfinite(Yiso))
    error('ntcp_isoE: el DVH contiene NaN/Inf.');
end
if ~(isfinite(TD50) && TD50 > 0)
    error('ntcp_isoE: TD50 debe ser finito y positivo.');
end
if ~(isfinite(m) && m > 0)
    error('ntcp_isoE: m debe ser finita y positiva.');
end
if ~(isfinite(a) && a >= 0)
    error('ntcp_isoE: a debe ser finito y no negativo.');
end
if a == 0 && ~(isfinite(ab) && ab > 0)
    error('ntcp_isoE: ab debe ser finito y positivo cuando a = 0.');
end

% DVH acumulativo [%] -> diferencial (fraccion de volumen por bin).
V = Yiso(:) / 100;
dV = -diff(V);
dV = [dV; 0];
dV = max(dV, 0);
if sum(dV) == 0
    dV = ones(size(Xiso)) / numel(Xiso);
else
    dV = dV / sum(dV);
end

if a == 0
    % Endpoint tipo Dmax (mucosa): prescripcion al maximo + EQD2 2 Gy/fx.
    Dmax = max(Xiso);
    Dref = Dmax * (ab + Dmax) / (ab + 2);
else
    % Endpoint tipo EUD (cerebro a=10): media generalizada sobre el DVH.
    Dref = (sum(dV .* (Xiso .^ a))) ^ (1 / a);
end

t = (Dref - TD50) / (m * TD50);
NTCP = 0.5 * (1 + erf(t / sqrt(2)));

if ~isfinite(NTCP) || NTCP < 0 || NTCP > 1
    error('ntcp_isoE: NTCP fuera de [0,1] o no finito (t=%.3g).', t);
end
end
