function T = rigidTransform(p)
%RIGIDTRANSFORM Construye la transformacion rigida homogenea 4x4.
%   T = RIGIDTRANSFORM(P) con P = [x y z alpha beta gamma]:
%     - P(1:3) : traslacion en mm
%     - P(4:6) : angulos de Euler [alpha beta gamma] en GRADOS
%   La transformacion mapea coordenadas del MARCO DEL PACIENTE al MARCO DEL
%   HAZ:  r_haz = T(1:3,1:3) * r_paciente + T(1:3,4).
%   La rotacion usa la convencion de eul2rot (R = Rz*Ry*Rx, grados).
%
%   Entrada : P [1x6] o [6x1] (mm, grados).
%   Salida  : T [4x4] (mm).
%
%   Ver tambien: eul2rot.m, transformPatient.m

p = p(:).';                       % normaliza a fila
if numel(p) ~= 6
    error('BNCT:rigidTransform:entradaInvalida', ...
        'P debe tener 6 elementos [x y z alpha beta gamma]; recibio %d.', numel(p));
end

T = eye(4);
T(1:3, 4)   = p(1:3);             % traslacion (mm)
T(1:3,1:3)  = eul2rot(p(4:6));    % rotacion (grados -> rad dentro de eul2rot)
end
