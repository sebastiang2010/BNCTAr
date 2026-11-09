function gen_golden()
% GEN_GOLDEN — Valores de referencia MATLAB para validar el puerto Python UTCP.
% Ejecutar: gen_golden  (desde cualquier carpeta; el script hace cd a su dir).
% Salida:   <here>/golden_utcp.json
%
% Usa: isodose.m y TCPc.m de la RAIZ del repo (solo lectura) y ntcp_isoE.m
% copiado de UTCP/utils a esta carpeta (UTCP original intacto).

here = fileparts(mfilename('fullpath'));           % .../SlicerModule/test/matlab_ref
repo = fileparts(fileparts(fileparts(here)));      % raiz del repo
addpath(repo);                                     % isodose.m, TCPc.m
addpath(here);                                     % ntcp_isoE.m (copia)

% ------------------------------------------------------------------
% Datos sinteticos deterministas (3x4x5), column-major al exportar
% ------------------------------------------------------------------
rng(42, 'twister');
sz = [3 4 5];
n = prod(sz);
raw.boro = 0.010 * rand(sz);
raw.thn  = 0.005 * rand(sz);
raw.fast = 0.008 * rand(sz);
raw.g    = 0.002 * rand(sz);
raw.boro(1,1,1) = 0; raw.thn(1,1,1) = 0;      % voxel sin dosis -> NaN
raw.fast(1,1,1) = 0; raw.g(1,1,1)   = 0;

b10 = 20; bratio = 3.5;
p14 = [0.35, 0.035, 3.090888558919538, 7.034e-17, ...
       2.468388470867692, 1.0205e-16, 2.468388470867692, 1.0205e-16, ...
       38.961, 216.4502, 0.8, 0.2, 0.2, 0.8];
tcp_params = [0.34, 8.3e-9, 22, 0.33, 34.632, 1212.1, 0.53, 0.47, 0.2, 0.8];
ntcp_mucosa = struct('TD50', 39.8, 'm', 0.17, 'a', 0, 'ab', 10);
ntcp_brain  = struct('TD50', 13.9049, 'm', 0.1766, 'a', 10, 'ab', 3);
tirs = [30, 45, 60];

% ------------------------------------------------------------------
% isodose.m por tiempo + TCPc.m + ntcp_isoE.m sobre el mismo DVH
% ------------------------------------------------------------------
D = struct();
D.boro = raw.boro * b10;      % isodose espera boro ya x ppm
D.thn = raw.thn; D.fast = raw.fast; D.g = raw.g;

dv_voxel = 0.027;             % cm^3
out = struct('dims', sz, 'b10', b10, 'bratio', bratio, ...
    'p14', {p14}, 'tcp_params', {tcp_params}, ...
    'ntcp_mucosa', {ntcp_mucosa}, 'ntcp_brain_ref', {ntcp_brain}, ...
    'dv_voxel', dv_voxel, 'tir', {[]}, 'diso', {{}}, ...
    'tcp', {[]}, 'peud', {[]}, 'ntcp_skin', {[]}, 'ntcp_brain', {[]}, ...
    'comps_raw', {{}});

out.comps_raw = {raw.boro(:), raw.thn(:), raw.fast(:), raw.g(:)};

for k = 1:numel(tirs)
    t = tirs(k);
    Diso = isodose(D, t, bratio, p14);
    % TCP sobre la muestra de voxels finitos (dv escalar, TCPc.m)
    Deq = Diso(~isnan(Diso));
    [TCP, PEUD] = TCPc(Deq, dv_voxel, tcp_params);
    % DVH acumulado del mismo conjunto (orden ascendente, Y % acumulado)
    Xiso = sort(Deq(:));
    N = numel(Xiso);
    Yiso = 100 * (1 - (0:N-1)' / N);
    n_skin = ntcp_isoE(Xiso, Yiso, ntcp_mucosa.TD50, ntcp_mucosa.m, ...
                       ntcp_mucosa.a, ntcp_mucosa.ab);
    n_brain = ntcp_isoE(Xiso, Yiso, ntcp_brain.TD50, ntcp_brain.m, ...
                        ntcp_brain.a, ntcp_brain.ab);
    out.tir{k} = t;
    out.diso{k} = Diso(:);     % column-major
    out.tcp{k} = TCP;
    out.peud{k} = PEUD;
    out.ntcp_skin{k} = n_skin;
    out.ntcp_brain{k} = n_brain;
end

txt = jsonencode(out, 'PrettyPrint', true);
fid = fopen(fullfile(here, 'golden_utcp.json'), 'w');
fwrite(fid, txt, 'char');
fclose(fid);
fprintf('golden_utcp.json escrito (%d tiempos)\n', numel(tirs));
end
