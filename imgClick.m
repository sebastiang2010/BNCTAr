function imgClick(hObject, eventdata)
global patient;
puntos = eventdata.IntersectionPoint;
patient.px = fix(puntos(2));
patient.py = fix(puntos(1));
plot_choosen_dose(hObject,0);
end

