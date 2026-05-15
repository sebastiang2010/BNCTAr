% Leer las imágenes de dosis y CT
dosis=patient.Dt.boro; 
ct=double(patient.stack); 


% Normalizar las imágenes en el rango [0, 1]
dosis = (dosis - min(dosis(:))) / (max(dosis(:)) - min(dosis(:)));
ct = (ct - min(ct(:))) / (max(ct(:)) - min(ct(:)));

% Especificar el optimizador y la métrica de similitud
optimizer = registration.optimizer.RegularStepGradientDescent;
metric = registration.metric.MattesMutualInformation;

% Crear la estructura de opciones
options = registration.optimizer.RegularStepGradientDescent();
options.MaximumIterations = 100;
options.MinimumStepLength = 1e-4;
options.MaximumStepLength = 0.0625;

% Aplicar la transformación afin para alinear las imágenes
tform = imregtform(dosis, ct, 'affine', optimizer, metric, 'PyramidLevels', 3, 'InitialTransformation', [], 'DisplayOptimization', []);

% Verificar que se ha obtenido una transformación válida
if isfield(tform,'T')
    dosis_aligned = imwarp(dosis, tform);
    % Asignar pesos a las imágenes antes de fusionarlas
    peso_dosis = 0.5;
    peso_ct = 0.5;

    % Fusionar las imágenes con ponderación de intensidad
    fusion = peso_dosis * dosis_aligned + peso_ct * ct;

    % Mostrar la imagen resultante
    imshow(fusion);
else
    disp('No se pudo obtener una transformación válida');
end
