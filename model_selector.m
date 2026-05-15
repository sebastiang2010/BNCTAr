% options:
% 1 - iso model
% 2 - tcp model
% 3 - ntcp model 
function parametros = model_selector(model_option)

fid = fopen('model_list.txt');
line = fgetl(fid);
while ischar(line)
    % salteo lineas comentadas y lineas en blanco
    if ((isempty(line)) || (length(line)>=1) && (line(1)=='%'))
        line = fgetl(fid);
        continue
    end
    
    if contains(line,'# ISOEFFECTIVE MODELS:')
        model_type = 1;
        line = fgetl(fid);
        continue
    end
    
    if contains(line,'# TCP MODELS:')
        model_type = 2;
        line = fgetl(fid);
        continue
    end
    
    if contains(line,'# NTCP MODELS:')
        model_type = 3;
        line = fgetl(fid);
        continue
    end

    if contains(line,'- MODEL')
        model_number = sscanf(line,'- MODEL %d:');
        line = fgetl(fid);
        model_name = line;
        
        if model_type == 1
            iso_model{model_number,1} = model_name;
        elseif model_type == 2
            tcp_model{model_number,1} = model_name;
        elseif model_type == 3
            ntcp_model{model_number,1} = model_name;
        end
        i=2;
        line = fgetl(fid);
        continue
    end
    
    B = regexp(line, '(?<name>\w+)\s*=\s*(?<value>[^,]+)', 'names');
    val = str2double(B.value);
%     if length(B)==1
%         val = str2double([B{1}]);
%     elseif length(B)==2
%         val = str2double([B{1} '.' B{2} ]);
%     elseif length(B)==3
%         val = str2double([B{1} '.' B{2} 'E' B{3}]);
%     end
    
    if model_type == 1
        iso_model{model_number,i} = val;
    elseif model_type == 2
        tcp_model{model_number,i} = val;
    elseif model_type == 3
        ntcp_model{model_number,i} = val;
    end
    
    i=i+1;
 
    % line
    line = fgetl(fid);
    
end

switch model_option
    case 1
        model = iso_model;
    case 2
        model = tcp_model;
    case 3
        model = ntcp_model;
    case 4
        prescrip_model = [{'Mucosa'} iso_model(1,2:15) ntcp_model(1,2:5)
            {'Brain'} iso_model(2,2:15) ntcp_model(2,2:5)];
        model = prescrip_model;
end

seleccion = DropMenu(model);


switch model_option
    case 1
        alpha_r = iso_model{seleccion,2};
        beta_r = iso_model{seleccion,3};
        alpha_boro = iso_model{seleccion,4};
        beta_boro = iso_model{seleccion,5};
        alpha_thn = iso_model{seleccion,6};
        beta_thn = iso_model{seleccion,7};
        alpha_fast = iso_model{seleccion,8};
        beta_fast = iso_model{seleccion,9};
        tof = iso_model{seleccion,10};
        tos = iso_model{seleccion,11};
        pf_g = iso_model{seleccion,12};
        ps_g = iso_model{seleccion,13};
        pf_bnct = iso_model{seleccion,14};
        ps_bnct = iso_model{seleccion,15};
        
        parametros = [alpha_r, beta_r, ...
                      alpha_boro, beta_boro, ...
                      alpha_thn, beta_thn, ...
                      alpha_fast, beta_fast, ...
                      tof, tos, pf_g, ps_g, pf_bnct, ps_bnct];
    case 2
        name = tcp_model{seleccion,1};
        Ar = tcp_model{seleccion,2};
        Br = tcp_model{seleccion,3};
        c1 = tcp_model{seleccion,4};
        c2 = tcp_model{seleccion,5};
        tof = tcp_model{seleccion,6};
        tos = tcp_model{seleccion,7};
        pf_g = tcp_model{seleccion,8};
        ps_g = tcp_model{seleccion,9};
        pf_bnct = tcp_model{seleccion,10};
        ps_bnct = tcp_model{seleccion,11};
        
        parametros = {name, Ar, Br, c1, c2, tof, tos, pf_g, ps_g, pf_bnct, ps_bnct };
    case 3
        name = ntcp_model{seleccion,1};
        TD50g = ntcp_model{seleccion,2};
        mg = ntcp_model{seleccion,3};
        ab = ntcp_model{seleccion,4};
        a = ntcp_model{seleccion,5};
        
        parametros = {name, TD50g, mg, ab, a};

    case 4
        alpha_r = prescrip_model{seleccion,2};
        beta_r = prescrip_model{seleccion,3};
        alpha_boro =prescrip_model{seleccion,4};
        beta_boro = prescrip_model{seleccion,5};
        alpha_thn = prescrip_model{seleccion,6};
        beta_thn = prescrip_model{seleccion,7};
        alpha_fast = prescrip_model{seleccion,8};
        beta_fast =  prescrip_model{seleccion,9};
        tof =  prescrip_model{seleccion,10};
        tos =  prescrip_model{seleccion,11};
        pf_g = prescrip_model{seleccion,12};
        ps_g = prescrip_model{seleccion,13};
        pf_bnct = prescrip_model{seleccion,14};
        ps_bnct = prescrip_model{seleccion,15};
        TD50g = prescrip_model{seleccion,16};
        mg = prescrip_model{seleccion,17};
        ab = prescrip_model{seleccion,18};
        a = prescrip_model{seleccion,19};

        parametros = [alpha_r, beta_r, ...
                      alpha_boro, beta_boro, ...
                      alpha_thn, beta_thn, ...
                      alpha_fast, beta_fast, ...
                      tof, tos, pf_g, ps_g, pf_bnct, ps_bnct,...
                      TD50g, mg, ab, a];        
end

clear B fid i iso_model line model_name model_number model_type selection tcp_model val ntcp_model prescrip_model
end

function seleccion = DropMenu(model)

f = DropMenu_create_GUI(model);
uiwait(f)
h = guidata(f);

seleccion = h.choices;
close(f)
end

function f = DropMenu_create_GUI(model)
ss = get(0,'ScreenSize');
ww = 250;
hh = 150;
xpp = (ss(3)/2)-(ww/2);
ypp = (ss(4)/2)-(hh/2);

f = figure('menu','none','Toolbar','none','Position', [xpp ypp ww hh]);
h = struct('f',f);

if (length(model)==15)
    model_txt = 'Select Isoeffective Dose model:';
elseif(length(model)==11)
    model_txt = 'Select TCP model:';
elseif(length(model)==5)
    model_txt = 'Select NTCP model:';
elseif(length(model)==19)
    model_txt = 'Select prescription tissue:';
end

h.txt_strig = uicontrol('Parent',f,...
    'Style','Text',...
    'Units','Normalized',...
    'Position',[0.1 0.65 0.8 0.2],...
    'String',model_txt);
h.button_ok = uicontrol('Parent',f,...
    'Style','pushbutton',...
    'Units','Normalized',...
    'Position',[0.1 0.1 0.3 0.20],...
    'String','OK',...
    'Callback',@DropMenu_Callback);
h.button_cancel = uicontrol('Parent',f,...
    'Style','pushbutton',...
    'Units','Normalized',...
    'Position',[0.6 0.1 0.3 0.20],...
    'String','cancel',...
    'Callback',@DropMenu_Callback);
h.list1 = uicontrol('Parent',f,...
    'Style','popupmenu',...
    'Units','Normalized',...
    'Position',[0.1 0.45 0.8 0.2],...
    'String',model(:,1));

% set(f,'CloseRequestFcn', @DropMenu_CloseFunction)
guidata(h.f,h)
end

function DropMenu_Callback(hObject,eventdata)
h = guidata(hObject);
switch hObject
    case h.button_ok
        h.choices = get(h.list1,'Value');
    case h.button_cancel
        h.choices = 0;
end
guidata(h.f,h)
uiresume(h.f)
end

function DropMenu_CloseFunction(hObject,eventdata)
h = guidata(hObject);
h.choices = [0];
guidata(h.f,h)
uiresume(h.f)
% delete(gcf)
end