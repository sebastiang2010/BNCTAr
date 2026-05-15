function [D,Dosis1,Dosis2] = f_leer_multicell_v2(sthn,sfast,sg)% Leo la mesh de MCNP.

% Selección de cantidad de campos
choice = listdlg('PromptString','Seleccione cantidad de campos:',...
                 'SelectionMode','single',...
                 'ListString',{'1 campo','2 campos'});

if isempty(choice)
    error('No se seleccionó ninguna opción');
end

% ================= CASO 1 CAMPO =================
if choice == 1
    
    Dosis1 = lecturaDosis;
    
    tir1 = 1.0;
    tir2 = 0.0;
    
    % Dosis2 vacía compatible
    Dosis2 = Dosis1;
    Dosis2.Matrix = zeros(size(Dosis1.Matrix));
    Dosis2.Error  = zeros(size(Dosis1.Error));

% ================= CASO 2 CAMPOS =================
else
    
    Dosis1 = lecturaDosis;
    
    prompt = {'Field 1 weight/time:'};
    dlg_title = 'Field 1 - Weight:';
    answ = inputdlg(prompt,dlg_title,1,{'30'});
    tir1 = str2double(answ{1});
    
    Dosis2 = lecturaDosis;
    
    prompt = {'Field 2 weight/time:'};
    dlg_title = 'Field 2 - Weight:';
    answ = inputdlg(prompt,dlg_title,1,{'30'});
    tir2 = str2double(answ{1});
    
end

% ================= PESOS =================
wgt1 = tir1/(tir1+tir2);
wgt2 = tir2/(tir1+tir2);


% ================= COMBINACIÓN =================
D1.boro(:,:,:)  = wgt1*Dosis1.Matrix(:,:,:,1) + wgt2*Dosis2.Matrix(:,:,:,1);
D1.fastn(:,:,:) = wgt1*Dosis1.Matrix(:,:,:,2) + wgt2*Dosis2.Matrix(:,:,:,2);
D1.thn(:,:,:)   = wgt1*Dosis1.Matrix(:,:,:,3) + wgt2*Dosis2.Matrix(:,:,:,3);
D1.g(:,:,:)     = wgt1*Dosis1.Matrix(:,:,:,4) + wgt2*Dosis2.Matrix(:,:,:,4);

% ================= EJES =================
D.xc = Dosis1.X0:Dosis1.StepX:Dosis1.Xf;
D.yc = Dosis1.Y0:Dosis1.StepY:Dosis1.Yf;
D.zc = Dosis1.Z0:Dosis1.StepZ:Dosis1.Zf;

% ================= ESCALADO =================
D.boro(:,:,:) = sthn  * D1.boro(:,:,:);
D.fast(:,:,:) = sfast * D1.fastn(:,:,:);
D.thn(:,:,:)  = sthn  * D1.thn(:,:,:);
D.g(:,:,:)    = sg    * D1.g(:,:,:);

% ================= ERROR =================
D.E(:,:,:,:) = ((Dosis1.Matrix(:,:,:,:) .* Dosis1.Error(:,:,:,:)) + ...
                (Dosis2.Matrix(:,:,:,:) .* Dosis2.Error(:,:,:,:))) ./ ...
               (Dosis1.Matrix(:,:,:,:) + Dosis2.Matrix(:,:,:,:));

end