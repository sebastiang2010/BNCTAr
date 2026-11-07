function R = eul2rot(angles)
%EUL2ROT Matriz de rotación 3x3 a partir de ángulos de Euler en grados.
%   R = EUL2ROT(ANGLES) con ANGLES = [alpha beta gamma] en GRADOS devuelve
%   la matriz de rotación activa R = Rz(gamma) * Ry(beta) * Rx(alpha)
%   (convención XYZ: primero rotación en X, luego Y, luego Z).
%
%   Entrada : ANGLES [1x3] grados.
%   Salida  : R [3x3] adimensional, ortonormal.
%
%   Ver tambien: rigidTransform.m

a = deg2rad(angles(1));   % rotacion en X (rad)
b = deg2rad(angles(2));   % rotacion en Y (rad)
c = deg2rad(angles(3));   % rotacion en Z (rad)

Rx = [1 0      0;
      0 cos(a) -sin(a);
      0 sin(a)  cos(a)];

Ry = [ cos(b) 0 sin(b);
       0      1 0;
      -sin(b) 0 cos(b)];

Rz = [cos(c) -sin(c) 0;
      sin(c)  cos(c) 0;
      0       0      1];

% Combinacion Z*Y*X (convencion XYZ del proyecto, NO modificar sin
% actualizar rigidTransform y la documentacion de marcos del README)
R = Rz * Ry * Rx;
end
