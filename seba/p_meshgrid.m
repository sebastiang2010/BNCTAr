% Definir vectores x e y
x = -3:0.5:3;
y = -2:0.5:2;

% Crear matrices X e Y
[X, Y] = meshgrid(x, y);

% Graficar la malla
figure
plot(X, Y, '.', 'MarkerSize', 10)