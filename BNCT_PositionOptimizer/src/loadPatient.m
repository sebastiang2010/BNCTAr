function patientStruct = loadPatient(maskPath)
%LOADPATIENT Carga una mascara NIfTI y arma la estructura del paciente.
%   patientStruct = LOADPATIENT(MASKPATH) lee un volumen NIfTI (.nii/.nii.gz),
%   construye la mascara binaria (valor > 0), extrae la superficie con
%   isosurface y devuelve:
%     .rawData    : volumen crudo tal como fue leido
%     .mask       : mascara binaria DOUBLE (0/1), tamano [nx ny nz]
%     .surface    : struct con .points (Nx3, mm, marco paciente) y .faces (Mx3)
%     .voxelSize  : [1x3] mm (de niftiinfo si esta disponible; default [1 1 1])
%     .origin     : [1x3] mm (offset de niftiinfo si esta disponible; default [0 0 0])
%     .gridX/.gridY/.gridZ : ejes ndgrid en mm (origin + (0:n-1)*voxelSize)
%
%   CONVENCION DE EJES: mask(i,j,k) con i<->x, j<->y, k<->z (formato ndgrid),
%   consistente con sampleDoseMaps, transformPatient y los mapas de dosis
%   PHITS. La superficie se devuelve en milimetros (marco del paciente).
%
%   PARA ARMAR UN PACIENTE COMPLETO (cuerpo + tumor + OARs), combine varias
%   llamadas (ver README):
%       p  = loadPatient('cuerpo.nii');
%       t  = loadPatient('tumor.nii');   p.tumorMask = t.mask;
%       m  = loadPatient('mucosa.nii');  p.oars.mucosa = struct('mask', m.mask);
%   Si una ROI tiene grilla propia distinta, guarde tambien sus .gridX/Y/Z
%   dentro de la entrada correspondiente de p.oars (o p.tumorGrid).
%
%   Requiere Image Processing Toolbox (niftiread/niftiinfo; niftiinfo es
%   opcional: si falla, se usan defaults y se emite una advertencia).

if exist(maskPath, 'file') ~= 2
    error('BNCT:loadPatient:archivoNoEncontrado', ...
        'No se encontro el archivo de mascara: %s', maskPath);
end

% ---------- lectura robusta del volumen ----------
try
    rawData = niftiread(maskPath);            % Image Processing Toolbox
catch errLectura
    error('BNCT:loadPatient:falloNiftiread', ...
        ['No se pudo leer "%s" con niftiread (Image Processing Toolbox).\n' ...
         'Verifique que el toolbox este instalado y que el archivo sea un NIfTI valido.\n' ...
         'Detalle original: %s'], maskPath, errLectura.message);
end

% ---------- mascara binaria (DOUBLE, no logical: isosurface/aritmetica) ----------
mask = double(rawData > 0);
if ~any(mask(:))
    warning('BNCT:loadPatient:mascaraVacia', ...
        'La mascara "%s" no tiene voxels con valor > 0.', maskPath);
end

% ---------- voxelSize y origin desde el encabezado NIfTI (best-effort) ----------
voxelSize = [1 1 1];      % mm
origin    = [0 0 0];      % mm
try
    info = niftiinfo(maskPath);
    md = info.Metadata;
    if isfield(md, 'PixDim') && numel(md.PixDim) >= 4
        sp = md.PixDim(2:4);                         % pixdim(2:4) = x,y,z (mm)
        if all(isfinite(sp)) && all(sp > 0)
            voxelSize = double(sp(:).');
        end
    end
    cand = [];
    if isfield(md, 'QoffsetX')
        cand = [md.QoffsetX, md.QoffsetY, md.QoffsetZ];
    elseif isfield(md, 'SrowX') && isfield(md, 'SrowY') && isfield(md, 'SrowZ')
        cand = [md.SrowX(4), md.SrowY(4), md.SrowZ(4)];
    end
    if ~isempty(cand) && all(isfinite(cand))
        origin = double(cand(1:3));
    end
catch
    warning('BNCT:loadPatient:sinNiftiinfo', ...
        ['niftiinfo no devolvio metadata util para "%s"; ' ...
         'se usan voxelSize=[1 1 1] mm y origin=[0 0 0].'], maskPath);
end

% ---------- superficie (isosurface opera en coordenadas de indice 1-based) ----------
% ATENCION: la forma de 2 salidas devuelve [CARAS, VERTICES] en ese orden.
[faces, verts] = isosurface(mask, 0.5);
verts = double(verts);
% conversion indice -> mm: x = origin + (i-1)*voxelSize
points = (verts - 1) .* voxelSize(:).' + origin(:).';

% ---------- estructura de salida ----------
patientStruct = struct();
patientStruct.rawData   = rawData;
patientStruct.mask      = mask;
patientStruct.surface   = struct('points', points, 'faces', double(faces));
patientStruct.voxelSize = voxelSize;
patientStruct.origin    = origin;
patientStruct.gridX     = origin(1) + (0:size(mask,1)-1) * voxelSize(1);
patientStruct.gridY     = origin(2) + (0:size(mask,2)-1) * voxelSize(2);
patientStruct.gridZ     = origin(3) + (0:size(mask,3)-1) * voxelSize(3);
end
