function  tr=f_prescripcion(h) 

prescripcion_porc=98; 
prescripcion_dosis=20; 

 prompt = {'Thermal:','Fast:','Gamma:'};
 dlg_title = 'Fraction 1 - Scaling Factors:';
 num_lines = 1;
 defaultans = {'0.97','1.05','1.07'};
 answ = inputdlg(prompt,dlg_title,num_lines,defaultans);


%h=gco;
porcentaje=h.YData; 
dosis=h.XData; 


porcen_new=porcentaje; 
porcen_new(end+1)=prescripcion_porc; 
porc_new=sort(porcen_new,'descend'); 

dosis2=interp1(porcentaje,dosis,porc_new);

ind=porc_new==prescripcion_porc; 

tasa_dosis=dosis2(ind); 

tr=prescripcion_dosis/tasa_dosis; 

end %funcion 