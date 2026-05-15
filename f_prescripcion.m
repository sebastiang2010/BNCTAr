function  tr=f_prescripcion(h) 

%prescripcion_porc=98;
%prescripcion_dosis=20;

prompt = {'Procentaje', 'Prescrpcion (Gy)'};
dlg_title = '';
num_lines = 1;
defaultans = {'98','20'};
answ = inputdlg(prompt,dlg_title,num_lines,defaultans);

prescripcion_porc=str2double(answ{1,1});
prescripcion_dosis=str2double(answ{2,1});

porcentaje=h.YData; 
dosis=h.XData; 

porcen_new=porcentaje; 
porcen_new(end+1)=prescripcion_porc; 
porcen_new=unique(porcen_new, 'stable');
porc_new=sort(porcen_new,'descend'); 

dosis2=interp1(porcentaje,dosis,porc_new);

ind=porc_new==prescripcion_porc; 

tasa_dosis=dosis2(ind); 

tr=prescripcion_dosis/tasa_dosis; 

end %funcion 