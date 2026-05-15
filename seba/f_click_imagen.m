function f_click_imagen(h) 
global p
puntos = eventdata.IntersectionPoint;
p.x = fix(puntos(2));
p.y = fix(puntos(1));
end