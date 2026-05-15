clc 

% Definir la geometría voxelizada
%A = randi([0 1], [20 20 20]); % Ejemplo de geometría voxelizada aleatoria
A=zeros(20,20,20); 
A(5:10,7:12,5:15)=1;

s=size(A); 

figure(1)
imshow(A(:,:,10),[])

% Crear un vector de coordenadas para cada voxel
[x,y,z] = meshgrid(1:s(2), 1:s(1), 1:s(3));

% Crear una matriz que contenga las coordenadas de todos los voxels
V = reshape([x(:), y(:), z(:)], [], 3);

% Encontrar los índices de los voxels que corresponden a la región de interés
I = find(A(:) == 1);

% Crear la triangulación de Delaunay utilizando solo los puntos correspondientes a la región de interés
DT = delaunayTriangulation(V(I,:));
T = DT.ConnectivityList;

% Visualizar la malla tetraédrica
%figure(2)
tetramesh(T, V);

%% 
% Definir las dimensiones de la geometría
x_range = [-5, 5];
y_range = [-7, 7];
z_range = [-10, 10];

% Definir la resolución de la geometría
dx = 0.25;
dy = 0.25;
dz = 0.25;

% Crear un vector de coordenadas para cada voxel
[x, y, z] = meshgrid(x_range(1):dx:x_range(2), y_range(1):dy:y_range(2), z_range(1):dz:z_range(2));

% Crear una matriz que contenga las coordenadas de todos los voxels
V = reshape([x(:), y(:), z(:)], [], 3);

% Crear la geometría de prueba
A = zeros(size(x));
A(5/dx:10/dx, 7/dy:12/dy, 5/dz:15/dz) = 1;

% Encontrar los índices de los voxels que corresponden a la región de interés
ind = find(A(:) == 1);

V1=V(ind,:);

% Crear la triangulación de Delaunay utilizando solo los puntos correspondientes a la región de interés
DT = delaunayTriangulation(V1);
T = DT.ConnectivityList;

% Visualizar la malla tetraédrica
figure(3)
tetramesh(T, V1);

figure(4)
tetramesh(DT)

