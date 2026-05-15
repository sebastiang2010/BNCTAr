function f_imgClick(hObject,eventdata)
     %puntos = eventdata.IntersectionPoint;
     punto = get(gca,'CurrentPoint');
     x = round(punto(1,1));  
     y = round(punto(1,2));     
     assignin('base', 'x', x);
     assignin('base', 'y', y);
     fig_handle = gcf;
     fig_num = findall(0,'Type','figure','NumberTitle','on');
     num_fig = find(fig_num == fig_handle);
     assignin('base','num_fig',num_fig); 
     
     drawnow
end