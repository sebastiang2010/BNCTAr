function [D,Dosis1,Dosis2] = leer_multicell(sthn,sfast,sg)% Leo la mesh de MCNP.

Dosis1 = lecturaDosis;
prompt = {'Field 1 weight/time:'};
dlg_title = 'Field 1 - Weight:';
num_lines = 1;
defaultans = {'30'};
answ = inputdlg(prompt,dlg_title,num_lines,defaultans);
        
tir1 = str2double(answ{1});

% Leo la mesh de MCNP.
Dosis2 = lecturaDosis;
prompt = {'Field 2 weight/time:'};
dlg_title = 'Field 2 - Weight:';
num_lines = 1;
defaultans = {'30'};
answ = inputdlg(prompt,dlg_title,num_lines,defaultans);
        
tir2 = str2double(answ{1});

%pesos de campos.
wgt1 = tir1/(tir2+tir1);
wgt2 = tir2/(tir2+tir1);

% Separo y peso componentes.
D1.boro(:,:,:) = wgt1*Dosis1.Matrix(:,:,:,1)+wgt2*Dosis2.Matrix(:,:,:,1);
D1.fastn(:,:,:) = wgt1*Dosis1.Matrix(:,:,:,2)+wgt2*Dosis2.Matrix(:,:,:,2);
D1.thn(:,:,:) = wgt1*Dosis1.Matrix(:,:,:,3)+wgt2*Dosis2.Matrix(:,:,:,3);
D1.g(:,:,:) = wgt1*Dosis1.Matrix(:,:,:,4)+wgt2*Dosis2.Matrix(:,:,:,4);

D.xc = Dosis1.X0:Dosis1.StepX:Dosis1.Xf;
D.yc = Dosis1.Y0:Dosis1.StepY:Dosis1.Yf;
D.zc = Dosis1.Z0:Dosis1.StepZ:Dosis1.Zf;

D.boro(:,:,:) = sthn*D1.boro(:,:,:);
D.fast(:,:,:) = sfast*D1.fastn(:,:,:);
D.thn(:,:,:) = sthn*D1.thn(:,:,:);
D.g(:,:,:) = sg*D1.g(:,:,:);
% el error esta calculado solo para la componente boro 
%D.E(:,:,:) = ((Dosis1.Matrix(:,:,:,1).*Dosis1.Error(:,:,:,1))+(Dosis2.Matrix(:,:,:,1).*Dosis2.Error(:,:,:,1)))./(Dosis1.Matrix(:,:,:,1)+Dosis2.Matrix(:,:,:,1));
D.E(:,:,:,:) = ((Dosis1.Matrix(:,:,:,:).*Dosis1.Error(:,:,:,:))+(Dosis2.Matrix(:,:,:,:).*Dosis2.Error(:,:,:,:)))./(Dosis1.Matrix(:,:,:,:)+Dosis2.Matrix(:,:,:,:));



%clear Dosis1 Dosis2