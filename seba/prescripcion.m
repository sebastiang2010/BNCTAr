clc 

prescripcion_porc=98; 
prescripcion_dosis=20; 

h=gco;
porcentaje=h.YData; 
dosis=h.XData; 


porcen_new=porcentaje; 
porcen_new(end+1)=prescripcion_porc; 
porc_new=sort(porcen_new,'descend'); 

dosis2=interp1(porcentaje,dosis,porc_new);

ind=find(porc_new==prescripcion_porc); 

tasa_dosis=dosis2(ind); 

tr=prescripcion_dosis/tasa_dosis; 