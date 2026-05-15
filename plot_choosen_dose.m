function plot_choosen_dose(hObject,plot_images)

global patient;


handles = guidata(hObject);


set(handles.txt_nstack,'String',[sprintf('%03i',patient.nimg) '/' sprintf('%03i',patient.stacksize)]);
set(handles.inp_xcor,'String',num2str(patient.px));
set(handles.inp_ycor,'String',num2str(patient.py));
set(handles.inp_zcor,'String',num2str(patient.nimg));

if isempty(patient.px) 
    sI=size(patient.stack);
    patient.px=round(sI(1)/2);
    clear sI
end 
if isempty(patient.py) 
    sI=size(patient.stack);
    patient.py=round(sI(2)/2);
    clear sI
end 
if isempty(patient.nimg) 
    sI=size(patient.stack);
    patient.nimg=round(sI(3)/2);
    clear sI
end 

gray_value = patient.stack(patient.px,patient.py,patient.nimg);
roi_name = char(patient.roi_list(patient.roi_level==gray_value));

switch patient.show
    case 0 %CT.
        txt = ['ROI: ' roi_name ' X: ' num2str(patient.px) '  Y: ' num2str(patient.py)];
        set(handles.txt_path,'String',txt);
    case 1 %Df.
        txt = ['ROI: ' roi_name ' X: ' num2str(patient.px) '  Y: ' num2str(patient.py) '  Dose: ' num2str(patient.Df(patient.px,patient.py,patient.nimg),'%.3f') ' Gy [ Tir: ' num2str(patient.tir,'%0.2f') 'min ][ NTCP: ' num2str(patient.ntcp) '| Diso:' num2str(patient.Deq(patient.px,patient.py,patient.nimg),'%.2f') '| b:' num2str(patient.rb(patient.px,patient.py,patient.nimg),'%.3f') ' n:' num2str(patient.rf(patient.px,patient.py,patient.nimg)+patient.rt(patient.px,patient.py,patient.nimg),'%.3f') ' g:' num2str(patient.rg(patient.px,patient.py,patient.nimg),'%.3f') ' ]' ];
        plot_dose = patient.Df;
        set(handles.txt_path,'String',txt);
    case 2 %Dw.
        txt = ['ROI: ' roi_name ' X: ' num2str(patient.px) '  Y: ' num2str(patient.py) '  Dose: ' num2str(patient.Dw(patient.px,patient.py,patient.nimg),'%.3f') ' Gy (RBE) [ Tir: ' num2str(patient.tir,'%0.2f') 'min ]'];
        plot_dose = patient.Dw;
        set(handles.txt_path,'String',txt);         
    case 3 %Diso.
        txt = ['ROI: ' roi_name ' X: ' num2str(patient.px) '  Y: ' num2str(patient.py) '  Dose: ' num2str(patient.Diso(patient.px,patient.py,patient.nimg),'%.3f') ' Gy (isoE) [ Tir: ' num2str(patient.tir,'%0.2f') 'min ]'];
        plot_dose = patient.Diso;
        set(handles.txt_path,'String',txt);
    case 4 %Relative Error.
        plot_dose = patient.Dt.E(:,:,:,patient.nE);
        a=plot_dose(patient.px,patient.py,patient.nimg)*100;
        %txt = ['ROI: ' roi_name ' X: ' num2str(patient.px) '  Y: ' num2str(patient.py) '  Relative Error: ' num2str(patient.Dt.E(patient.px,patient.py,patient.nimg),'%.2f') ];
        txt = ['ROI: ' roi_name ' X: ' num2str(patient.px) '  Y: ' num2str(patient.py) '  Relative Error (%): ' num2str(a, '%.5f') ];
        set(handles.txt_path,'String',txt);
        txt=['Componet:  ' num2str(patient.nE)];
        set(handles.txt_status,'String',txt)
        clear a 
end

if(plot_images)
    % Plot CT Scan.
    I = im2double(patient.stack(:,:,patient.nimg));
    [m,n] = size(I);
    levels = hist(reshape(I,[1 m*n]),patient.stacksizeXY);
    levels = find(levels)-1;
    

% %no entiendo para que esta          
%     % Plot ROI contours.
%     for i = 2:length(levels)
%         BW = false(m,n);
%         BW(patient.stack(:,:,patient.nimg)==levels(i))=true;
%         B = bwboundaries(BW,'noholes');
%         for j = 1:length(B)
%             bound = B{j};
%             I(sub2ind([m n],bound(:,1),bound(:,2))) = 0;
%         end
%     end
    
    Irgb = cat(3,I,I,I);
    if(~isfield(patient,'himg'))
        patient.himg = image(Irgb);
        axis off
    else
        set(patient.himg,'CData',Irgb)
    end
    
end     
    
set(patient.himg,'ButtonDownFcn',@imgClick)
    
hold all
 


if  (patient.show)
    max_dose = max(max(max(plot_dose())));
    min_dose = min(min(min(plot_dose())));
    %pmax = ind2sub(size(plot_dose),find(plot_dose==max(max(max(plot_dose)))));
    plot_dose = plot_dose(:,:,patient.nimg);
    %Plot Dose.
    if(~isfield(patient,'dose_surf'))
        patient.dose_surf = surf(plot_dose,'FaceAlpha',0.4);
        set(patient.dose_surf,'ButtonDownFcn',@imgClick)
        shading interp;
        colormap jet
    else
        set(patient.dose_surf,'ZData',plot_dose)
    end
    clim([min_dose max_dose]);
    
 

    if patient.show==4 
        if(~isfield(patient,'hisoc') && patient.show)
            [C,patient.hisoc] = contour(100.*(plot_dose), 'LineColor', 'k');
            clabel(C,patient.hisoc,'Color','white')
            patient.hisoc.ContourZLevel = 200;
            set(patient.hisoc,'ButtonDownFcn',@imgClick)
        elseif(patient.show)
            patient.hisoc.ZData = plot_dose.*100;
        end
    else
        % Plot Isocurves.
        if(~isfield(patient,'hisoc') && patient.show)
            [C,patient.hisoc] = contour(100.*(plot_dose./max_dose), 'LineColor', 'k');
            clabel(C,patient.hisoc,'Color','white')
            patient.hisoc.ContourZLevel = 200;
            set(patient.hisoc,'ButtonDownFcn',@imgClick)
        elseif(patient.show)
            patient.hisoc.ZData = 100.*(plot_dose./max_dose);
        end
    end

    if(patient.show==4)
        txt2 = 'min';
    else
        txt2 = 'max';
    end
    
    if(~isfield(patient,'hplotmax'))
        patient.hplotmax = plot3(patient.ydmax,patient.xdmax,200,'k*');
        set(patient.hplotmax,'ButtonDownFcn',@imgClick)
        patient.hplotmax_text = text(patient.ydmax+3,patient.xdmax+3,200,txt2,'Color','k');
        set(patient.hplotmax_text,'ButtonDownFcn',@imgClick)
    else
        set(patient.hplotmax,'XData',patient.ydmax,'YData',patient.xdmax)
        set(patient.hplotmax_text,'String',txt2)
        set(patient.hplotmax_text,'Position',[patient.ydmax+3 patient.xdmax+3 200])
    end


    if(patient.zdmax == patient.nimg)
        set(patient.hplotmax,'Visible','on')
        set(patient.hplotmax_text,'Visible','on')
    else
        set(patient.hplotmax,'Visible','off')
        set(patient.hplotmax_text,'Visible','off')
    end

end 

%% agrego las lineas 
if (isfield(patient,'p1cursor'))
    delete(patient.p1cursor);
    delete(patient.p2cursor);
    delete(patient.cursor);
end 
hold on
patient.p1cursor = plot3([patient.py patient.py] ,[1 512],[200 200],'w');
patient.p2cursor = plot3([1 512],[patient.px patient.px],[200 200],'w');
patient.cursor = plot3(patient.py,patient.px,200,'wo');
patient.p1cursor.Color(4) = 0.3;
patient.p2cursor.Color(4) = 0.3;
%patient.text = text(5,patient.stacksizeXY-10,100,txt,'Color','w');
set(patient.p1cursor,'ButtonDownFcn',@imgClick)
set(patient.p2cursor,'ButtonDownFcn',@imgClick)
hold off 

%% cortes sagital y coronal 
  
 f_cor_sag(patient)
 
%% analizo y grabo 
% if patient.info==0 && patient.show>=1
%     %patient.componet_info=f_evaluacion(patient);  
%     directorio=f_creo_directorio; 
%     archivo=[directorio,'info_patient']; 
% 
%     [patient.info_error,patient.info_error_roi]=f_error(patient);
%     patient.info_mesh=f_analisis_mesh(patient);
%     save(archivo,'patient')
%     patient.info=1; 
%     clc
%     %disp('   ')
%     %disp('Se grabo el archivo info_patient en: ')
%     %disp(directorio)
%     txt=['The file info_patient was saved in:  ',directorio];
%     msgbox(txt);
%  end 

  

end  


