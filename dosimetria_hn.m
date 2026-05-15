function varargout = dosimetria_hn(varargin)
% DOSIMETRIA_HN MATLAB code for dosimetria_hn.fig
%      DOSIMETRIA_HN, by itself, creates a new DOSIMETRIA_HN or raises the existing
%      singleton*.
%
%      H = DOSIMETRIA_HN returns the handle to a new DOSIMETRIA_HN or the handle to
%      the existing singleton*.
%
%      DOSsfasgvasgasgdsfdgsdIMETRIA_HN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DOSIMETRIA_HN.M with the given input arguments.
%
%      DOSIMETRIA_HN('Property','Value',...) creates a new DOSIMETRIA_HN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dosimetria_hn_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dosimetria_hn_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dosimetria_hn

% Last Modified by GUIDE v2.5 26-Oct-2022 18:36:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @dosimetria_hn_OpeningFcn, ...
    'gui_OutputFcn',  @dosimetria_hn_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before dosimetria_hn is made visible.
function dosimetria_hn_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dosimetria_hn (see VARARGIN)

% Choose default command line output for dosimetria_hn
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes dosimetria_hn wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = dosimetria_hn_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% btn_new_case: ininizialize new patient case.
function btn_new_case_Callback(hObject, eventdata, handles)

global patient;
patient=struc([]);

clc
if ishandle(2);close(2);end
if ishandle(3);close(3);end

set(handles.inp_xmin,'Enable','off');
set(handles.inp_ymin,'Enable','off');
set(handles.inp_zmin,'Enable','off');
set(handles.inp_xmax,'Enable','off');
set(handles.inp_ymax,'Enable','off');
set(handles.inp_zmax,'Enable','off');
set(handles.inp_xstep,'Enable','off');
set(handles.inp_ystep,'Enable','off');
set(handles.inp_zstep,'Enable','off');
set(handles.tbl_roi,'Visible','off');
set(handles.btn_left_arrow,'Enable','off');
set(handles.btn_rigth_arrow,'Enable','off');
set(handles.txt_nstack,'Enable','off');
set(handles.inp_xcor,'Enable','off');
set(handles.inp_ycor,'Enable','off');
set(handles.inp_zcor,'Enable','off');
set(handles.btn_go,'Enable','off');
set(handles.btn_load_1fr,'Enable','off');
set(handles.btn_load_2fr,'Enable','off');
set(handles.inp_tir,'Enable','off')
set(handles.btn_showMCNPerror,'Enable','off');
set(handles.btn_abs_dose,'Enable','off');
set(handles.btn_w_dose,'Enable','off');
set(handles.btn_isoe_dose,'Enable','off');
%set(handles.pnl_matrixrotation,'Enable','off');
set(handles.pnl_matrixrotation,'Enable','off');
set(handles.btn_mX,'Enable','off');
set(handles.btn_mY,'Enable','off');
set(handles.btn_mZ,'Enable','off');
set(handles.btn_sXY,'Enable','off');
set(handles.btn_sYZ,'Enable','off');
set(handles.btn_sZX,'Enable','off');
set(handles.btn_create_dvh,'Enable','off');
set(handles.btn_utcp,'Enable','off');
set(handles.btn_tir_point,'Enable','off');
%set(handles.btn_tir_region,'Enable','off');
set(handles.inp_xmin,'String','');
set(handles.inp_ymin,'String','');
set(handles.inp_zmin,'String','');
set(handles.inp_xmax,'String','');
set(handles.inp_ymax,'String','');
set(handles.inp_zmax,'String','');
set(handles.inp_xstep,'String','');
set(handles.inp_ystep,'String','');
set(handles.inp_zstep,'String','');
set(handles.inp_1stb10,'String','');
set(handles.inp_tir,'String','1')
set(handles.inp_xcor,'String','');
set(handles.inp_ycor,'String','');
set(handles.inp_zcor,'String','');
set(handles.txt_nstack,'String','');
set(handles.rbtn_showMCNPerror,'Value',0)
set(handles.rbtn_showDiso,'Value',0)
set(handles.rbtn_showDf,'Value',0)
set(handles.rbtn_showDw,'Value',0)

ax=gca;
cla(ax)
ax.Toolbar.Visible = 'on';

patient.patient_id = inputdlg('Patient ID');
set(handles.txt_patient_id,'String',[patient.patient_id{1}]);
set(handles.txt_patient_id,'Visible','on')
set(handles.btn_save_case,'Enable','on');
set(handles.btn_load_images,'Enable','on');

if isempty(patient.patient_id{1})
    txt='Unidentified patient';
    set(handles.txt_patient_id,'String',txt);
    set(handles.txt_patient_id,'ForegroundColor',[1 0 0]) %rojo
end
patient.info=0;


% btn_load_case: load existing patient case.
function btn_load_case_Callback(hObject, eventdata, handles)
global patient;

clc
clear patient
if ishandle(2);close(2);end
if ishandle(3);close(3);end

set(handles.inp_xmin,'Enable','off');
set(handles.inp_ymin,'Enable','off');
set(handles.inp_zmin,'Enable','off');
set(handles.inp_xmax,'Enable','off');
set(handles.inp_ymax,'Enable','off');
set(handles.inp_zmax,'Enable','off');
set(handles.inp_xstep,'Enable','off');
set(handles.inp_ystep,'Enable','off');
set(handles.inp_zstep,'Enable','off');
set(handles.tbl_roi,'Visible','off');
set(handles.btn_left_arrow,'Enable','off');
set(handles.btn_rigth_arrow,'Enable','off');
set(handles.txt_nstack,'Enable','off');
set(handles.inp_xcor,'Enable','off');
set(handles.inp_ycor,'Enable','off');
set(handles.inp_zcor,'Enable','off');
set(handles.btn_go,'Enable','off');
set(handles.btn_load_1fr,'Enable','off');
set(handles.btn_load_2fr,'Enable','off');
set(handles.inp_tir,'Enable','off')
set(handles.btn_showMCNPerror,'Enable','off');
set(handles.btn_abs_dose,'Enable','off');
set(handles.btn_w_dose,'Enable','off');
set(handles.btn_isoe_dose,'Enable','off');
%set(handles.pnl_matrixrotation,'Enable','off');
set(handles.pnl_matrixrotation,'Enable','off');
set(handles.btn_mX,'Enable','on');
set(handles.btn_mY,'Enable','on');
set(handles.btn_mZ,'Enable','on');
set(handles.btn_sXY,'Enable','on');
set(handles.btn_sYZ,'Enable','on');
set(handles.btn_sZX,'Enable','on');
set(handles.btn_create_dvh,'Enable','off');
set(handles.btn_utcp,'Enable','off');
set(handles.btn_tir_point,'Enable','off');
%set(handles.btn_tir_region,'Enable','off');
set(handles.inp_xmin,'String','');
set(handles.inp_ymin,'String','');
set(handles.inp_zmin,'String','');
set(handles.inp_xmax,'String','');
set(handles.inp_ymax,'String','');
set(handles.inp_zmax,'String','');
set(handles.inp_xstep,'String','');
set(handles.inp_ystep,'String','');
set(handles.inp_zstep,'String','');
set(handles.inp_1stb10,'String','');
set(handles.inp_tir,'String','1')
set(handles.inp_xcor,'String','');
set(handles.inp_ycor,'String','');
set(handles.inp_zcor,'String','');
set(handles.txt_nstack,'String','');
set(handles.rbtn_showMCNPerror,'Value',0)
set(handles.rbtn_showDiso,'Value',0)
set(handles.rbtn_showDf,'Value',0)
set(handles.rbtn_showDw,'Value',0)

ax=gca;
cla(ax)
ax.Toolbar.Visible = 'on';



set(handles.txt_status,'String','Loading Case...');
[filename, pathname] = uigetfile('*.mat', 'Select case mat file');
if isequal(filename,0)||isequal(pathname,0)
else
    clear patient
    load([pathname filename]);
    set(handles.txt_path,'String',['Path: ', pathname filename]);
end
set(handles.txt_patient_id,'String',[patient.patient_id{1}]);
set(handles.txt_patient_id,'Visible','on')
set(handles.btn_save_case,'Enable','on');
set(handles.btn_load_images,'Enable','off');

%borro cosas que no necesito
if (isfield(patient,'himg'))
    patient = rmfield(patient,'himg');
end
%if (isfield(patient,'info'))
%    patient = rmfield(patient,'info');
%end
patient.info=0;


if (isfield(patient,'stack'))
    set(handles.img,'Visible','on');
    axes(handles.img);
    patient.nimg = fix(patient.stacksize/2);
    patient.show = 0;

    n_roi = length(patient.roi_level);
    data = cell(n_roi, 4);
    if ~isfield(patient,'roi_tumor')
        patient.roi_tumor = zeros(n_roi,1);
    end
    if ~isfield(patient,'roi_color')
        patient.roi_color = lines(n_roi);
    end
    cmap = patient.roi_color;
    for i=1:n_roi
        data{i,1} = ['<html><table border=0 width=300 bgcolor=' rgb2hex([patient.roi_level(i) patient.roi_level(i) patient.roi_level(i)]) '><TR><TD> </TD></TR></table></html>'];
        data{i,2} = char(patient.roi_list(i));
        data{i,3} = patient.dvh_roi(i);
        data{i,4} = patient.roi_tumor(i);
    end

    set(handles.tbl_roi,'data',data);
    set(handles.tbl_roi,'ColumnName',{'Color','ROI Name','DVH','Tumor'});
    set(handles.tbl_roi,'ColumnWidth',{30,100,30,40});
    set(handles.tbl_roi,'ColumnEditable',[false true true true]);
    set(handles.tbl_roi,'Visible','on');

    ax=gca;
    cla(ax)
    ax.Toolbar.Visible = 'on';

    plot_choosen_dose(hObject,1);

    set(handles.btn_left_arrow,'Enable','on');
    set(handles.btn_rigth_arrow,'Enable','on');
    set(handles.txt_nstack,'Enable','on');

    set(handles.inp_xcor,'Enable','on');
    set(handles.inp_ycor,'Enable','on');
    set(handles.inp_zcor,'Enable','on');
    set(handles.btn_go,'Enable','on');

    set(handles.btn_load_1fr,'Enable','on');
    set(handles.btn_load_2fr,'Enable','on');

    set(handles.inp_xmin,'Enable','on');
    set(handles.inp_ymin,'Enable','on');
    set(handles.inp_zmin,'Enable','on');
    set(handles.inp_xmax,'Enable','on');
    set(handles.inp_ymax,'Enable','on');
    set(handles.inp_zmax,'Enable','on');
    set(handles.inp_xstep,'Enable','on');
    set(handles.inp_ystep,'Enable','on');
    set(handles.inp_zstep,'Enable','on');


end
set(handles.btn_showMCNPerror,'Enable','on');
set(handles.btn_abs_dose,'Enable','on');
set(handles.btn_w_dose,'Enable','on');
set(handles.btn_isoe_dose,'Enable','on');
set(handles.btn_tir_point,'Enable','on');
set(handles.inp_tir,'Enable','on');
set(handles.inp_1stb10,'Enable','on')



if (isfield(patient,'boron1'))
    set(handles.inp_1stb10,'String',patient.boron1);
end
if (isfield(patient,'xmin'))
    set(handles.inp_xmin,'String',patient.xmin);
end
if (isfield(patient,'ymin'))
    set(handles.inp_ymin,'String',patient.ymin);
end
if (isfield(patient,'zmin'))
    set(handles.inp_zmin,'String',patient.zmin);
end
if (isfield(patient,'xmax'))
    set(handles.inp_xmax,'String',patient.xmax);
end
if (isfield(patient,'ymax'))
    set(handles.inp_ymax,'String',patient.ymax);
end
if (isfield(patient,'zmax'))
    set(handles.inp_zmax,'String',patient.zmax);
end
if (isfield(patient,'xstep'))
    set(handles.inp_xstep,'String',patient.xstep);
end
if (isfield(patient,'ystep'))
    set(handles.inp_ystep,'String',patient.ystep);
end
if (isfield(patient,'zstep'))
    set(handles.inp_zstep,'String',patient.zstep);
end


set(handles.txt_status,'String','Load Case...OK');

% btn_save_case: save current patient case.
function btn_save_case_Callback(hObject, eventdata, handles)
global patient;

[filename, pathname] = uiputfile('*.mat', 'Save case as...');
save([pathname filename],'patient');
%txt=['Path: ',pathname,filename];
%set(handles.txt_path,'String',txt);
msgbox(' Save OK', 'Confirmation', 'help');

% btn_load_images: load new patient images.
function btn_load_images_Callback(hObject, eventdata, handles)

global patient;

[file,dir] = uigetfile('*.tif;*.tiff','Select tiff file');
if isequal(file,0)||isequal(dir,0)
    patient.stack=[];
else
    fimg = fullfile(dir,file);
    patient.stacksize = size(imfinfo(fimg),1);
    GralInf=imfinfo(fimg);
    patient.stacksizeXY=GralInf.Width;
    SizeXY=patient.stacksizeXY;
    for i=1:patient.stacksize
        patient.stack(:,:,i)=imread(fimg,'tif',i);
    end
    set(handles.img,'Visible','on');
    axes(handles.img);
    patient.nimg = fix(patient.stacksize/2);
    patient.px = size(patient.stack,1);
    patient.py = size(patient.stack,2);
    levels = histcounts(reshape(double(patient.stack),1,numel(patient.stack)),-0.5:255.5);
    levels = find(levels)-1;

    n_roi = length(levels);
    data = cell(n_roi, 4);
    cmap = lines(n_roi);
    for i=1:n_roi
        data{i,1} = ['<html><table border=0 width=300 bgcolor=' rgb2hex([levels(i) levels(i) levels(i)]) '><TR><TD> </TD></TR></table></html>'];
        data{i,2} = i;
        data{i,3} = 0;
        data{i,4} = 0;
    end

    patient.roi_list = string(1:n_roi);
    patient.roi_level = levels;
    patient.dvh_roi = zeros(n_roi,1);
    patient.roi_tumor = zeros(n_roi,1);
    patient.roi_color = cmap;

    set(handles.tbl_roi,'data',data);
    set(handles.tbl_roi,'ColumnName',{'Color','ROI Name','DVH','Tumor'});
    set(handles.tbl_roi,'ColumnWidth',{30,100,30,40});
    set(handles.tbl_roi,'ColumnEditable',[false true true true]);
    set(handles.btn_left_arrow,'Enable','on');
    set(handles.btn_rigth_arrow,'Enable','on');

    set(handles.inp_xcor,'Enable','on');
    set(handles.inp_ycor,'Enable','on');
    set(handles.inp_zcor,'Enable','on');
    set(handles.btn_go,'Enable','on');

    set(handles.txt_nstack,'Enable','on');
    set(handles.btn_load_1fr,'Enable','on');

    set(handles.inp_xmin,'Enable','on');
    set(handles.inp_ymin,'Enable','on');
    set(handles.inp_zmin,'Enable','on');
    set(handles.inp_xmax,'Enable','on');
    set(handles.inp_ymax,'Enable','on');
    set(handles.inp_zmax,'Enable','on');
    set(handles.inp_xstep,'Enable','on');
    set(handles.inp_ystep,'Enable','on');
    set(handles.inp_zstep,'Enable','on');
    set(handles.tbl_roi,'Visible','on');



    patient.show = 0;
    plot_choosen_dose(hObject,1);
    ax=gca;
    ax.Toolbar.Visible = 'on';
end

% inp_xmin: xmin input.
function inp_xmin_Callback(hObject, eventdata, handles)

global patient;

patient.xmin = str2double(get(hObject,'String'));

function inp_xmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% inp_ymin: ymin input.
function inp_ymin_Callback(hObject, eventdata, handles)

global patient;

patient.ymin = str2double(get(hObject,'String'));


function inp_ymin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% inp_zmin: zmin input.
function inp_zmin_Callback(hObject, eventdata, handles)

global patient;

patient.zmin = str2double(get(hObject,'String'));

function inp_zmin_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% inp_xmax: xmax input.
function inp_xmax_Callback(hObject, eventdata, handles)

global patient;

patient.xmax = str2double(get(hObject,'String'));

function inp_xmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% inp_ymax: ymax input.
function inp_ymax_Callback(hObject, eventdata, handles)

global patient;

patient.ymax = str2double(get(hObject,'String'));

function inp_ymax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% inp_zmax: zmax input.
function inp_zmax_Callback(hObject, eventdata, handles)

global patient;

patient.zmax = str2double(get(hObject,'String'));

function inp_zmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function inp_zstep_Callback(hObject, eventdata, handles)

global patient;

patient.zstep = str2double(get(hObject,'String'));


function inp_zstep_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function inp_ystep_Callback(hObject, eventdata, handles)

global patient;

patient.ystep = str2double(get(hObject,'String'));

function inp_ystep_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function inp_xstep_Callback(hObject, eventdata, handles)

global patient;

patient.xstep = str2double(get(hObject,'String'));

function inp_xstep_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in btn_left_arrow.
function btn_left_arrow_Callback(hObject, eventdata, handles)

global patient;

if(patient.nimg>1)
    patient.nimg =patient.nimg-1;
    plot_choosen_dose(hObject,1);
end

% --- Executes on button press in btn_rigth_arrow.
function btn_rigth_arrow_Callback(hObject, eventdata, handles)

global patient;

if(patient.nimg<patient.stacksize)
    patient.nimg =patient.nimg+1;
    plot_choosen_dose(hObject,1);
end

% --- Executes on button press in btn_load_1fr.
function btn_load_1fr_Callback(hObject, eventdata, handles)

global patient;

set(handles.txt_status,'String','Loading Dose Matrix...')

choice = questdlg('Select type of Dose Matrix:',...
                  'Dose Matrix',...
                  'Sera','Multicell','PHITS','Sera');

switch choice

    % ==========================
    case 'Sera'
    % ==========================
        [filename, pathname] = uigetfile('*.dat','Select Sera Dose Matrix');
        
        if ~(isequal(filename,0) || isequal(pathname,0))
            
            prompt = {'Neutrons:','Gamma'};
            dlg_title = 'Fraction 1 - Scaling Factors:';
            answ = inputdlg(prompt,dlg_title,1);

            sf_boron = str2double(answ{1});
            sf_g = str2double(answ{2});

            patient.D1c = leer_sera(fullfile(pathname,filename),sf_g,sf_boron);
            patient.D1c.E = zeros(size(patient.D1c.boro));
        end


    % ==========================
    case 'Multicell'
    % ==========================
        prompt = {'Thermal:','Fast:','Gamma:'};
        dlg_title = 'Fraction 1 - Scaling Factors:';
        defaultans = {'1','1','1'};
        
        answ = inputdlg(prompt,dlg_title,1,defaultans);

        sthn  = str2double(answ{1});
        sfast = str2double(answ{2});
        sg    = str2double(answ{3});

        choice2 = questdlg('Select type of MCNP output:',...
                           'MCNP output',...
                           'rssa (one file)',...
                           'sdef (two files)',...
                           'Dose .mat',...
                           'rssa (one file)');

        switch choice2

            case 'rssa (one file)'
                [patient.D1c,patient.beam.n1,patient.beam.n2] = ...
                    f_leer_multicell_v2(sthn,sfast,sg);

            case 'sdef (two files)'
                patient.D1c = leer_multicell2(sthn,sfast,sg);

            case 'Dose .mat'
                tipoarchivo = '.mat';
                currentdirectory = pwd;

                [archivo,directorio] = uigetfile(tipoarchivo,'Dose matrix');

                if ~(isequal(archivo,0) || isequal(directorio,0))

                    cd(directorio);

                    file = [directorio,archivo];
                    a = load(file);

                    patient.D1c = a.D1c;
                    clear a

                    cd(currentdirectory);
                end
        end


    % ==========================
    case 'PHITS'
    % ==========================
        for i = 1:2
            if i == 1
                D_n = f_lectura_phits_neutron();
                patient.D1c=D_n;     
            else
                % lectura de gamma 
                D_g=f_lectura_phits_gamma(); 
                patient.D1c.g = permute(D_g.Dg, [3 1 2]).*0.0; % esto deberia estar en la funcion 
                
            end
        end
        
        clear D_n D_g 
    
end


prompt = {'Blood B10 Concentration:'};
dlg_title = 'Fraction 1';
num_lines = 1;
defaultans = {'15'};
answ = inputdlg(prompt,dlg_title,num_lines,defaultans);

patient.boron1 = str2double(answ{1});

set(handles.inp_1stb10,'Enable','on');
set(handles.inp_1stb10,'String',num2str(patient.boron1));

set(handles.btn_showMCNPerror,'Enable','on');
set(handles.btn_abs_dose,'Enable','on');
set(handles.btn_w_dose,'Enable','on');
set(handles.btn_isoe_dose,'Enable','on');

a=get(handles.inp_xmin);
xmin=str2double(a.String);
a=get(handles.inp_ymin);
ymin=str2double(a.String);
a=get(handles.inp_zmin);
zmin=str2double(a.String);

a=get(handles.inp_xmax);
xmax=str2double(a.String);
a=get(handles.inp_ymax);
ymax=str2double(a.String);
a=get(handles.inp_zmax);
zmax=str2double(a.String);



a=get(handles.inp_xstep);
xstep=str2double(a.String);
a=get(handles.inp_ystep);
ystep=str2double(a.String);
a=get(handles.inp_zstep);
zstep=str2double(a.String);

xmp = xstep/2;
ymp = ystep/2;
zmp = zstep/2;

x = xmin+xmp:xstep:xmax-xmp;
y = ymin+ymp:ystep:ymax-ymp;
z = zmin+zmp:zstep:zmax-zmp;

%x=patient.D1c.xc;
%y=patient.D1c.yc;
%z=patient.D1c.zc;

%% agregar un chequeo para verificar que son igules
[X,Y,Z] = ndgrid(x,y,z);
[Xc,Yc,Zc] = ndgrid(patient.D1c.xc,patient.D1c.yc,patient.D1c.zc);

% Interpolo componente a componente.
set(handles.txt_status,'String','Interpolating...')


%patient.D1c.boro=permute(patient.D1c.boro,[2,1,3]);

patient.D1.boro=[];
patient.D1.fast=[];
patient.D1.thn=[];
patient.D1.g=[];
patient.D1.E=[];




patient.D1.boro = interpn(Xc,Yc,Zc,patient.D1c.boro,X,Y,Z,'linear');
% patient.D1.fast = interpn(Xc,Yc,Zc,patient.D1c.fast.*0,X,Y,Z,'linear');
patient.D1.fast= interpn(Xc,Yc,Zc,patient.D1c.fast,X,Y,Z,'linear');
patient.D1.thn = interpn(Xc,Yc,Zc,patient.D1c.thn,X,Y,Z,'linear');
patient.D1.g = interpn(Xc,Yc,Zc,patient.D1c.g,X,Y,Z,'linear');
% for i=1:4
%    patient.D1.E(:,:,:,i) = interpn(Xc,Yc,Zc,patient.D1c.E(:,:,:,i),X,Y,Z,'linear');
% end


patient.Dt.boro=[];
patient.Dt.fast=[];
patient.Dt.thn=[];
patient.Dt.g=[];
patient.Dt.E=[];

patient.Dt.boro=patient.D1.boro;
patient.Dt.fast = patient.D1.fast;
patient.Dt.thn = patient.D1.thn;
patient.Dt.g = patient.D1.g;
patient.Dt.E=patient.D1.E;

patient.tir = 1;


patient.xmin=xmin;
patient.ymin=ymin;
patient.zmin=zmin;
patient.xmax=xmax;
patient.ymax=ymax;
patient.zmax=zmax;
patient.xstep=xstep;
patient.ystep=ystep;
patient.zstep=zstep;




set(handles.btn_mX,'Enable','on');
set(handles.btn_mY,'Enable','on');
set(handles.btn_mZ,'Enable','on');
set(handles.btn_sXY,'Enable','on');
set(handles.btn_sYZ,'Enable','on');
set(handles.btn_sZX,'Enable','on');
set(handles.btn_load_2fr,'Enable','off');
set(handles.inp_tir,'Enable','on')
set(handles.txt_status,'String','Dose Matrix...OK')

set(handles.inp_xmin,'Enable','off');
set(handles.inp_ymin,'Enable','off');
set(handles.inp_zmin,'Enable','off');
set(handles.inp_xmax,'Enable','off');
set(handles.inp_ymax,'Enable','off');
set(handles.inp_zmax,'Enable','off');
set(handles.inp_xstep,'Enable','off');
set(handles.inp_ystep,'Enable','off');
set(handles.inp_zstep,'Enable','off');
set(handles.btn_load_1fr,'Enable','off')
set(handles.btn_load_images,'Enable','off')
set(handles.inp_1stb10,'Enable','on')


size_stack = size(patient.stack);
size_boro  = size(patient.Dt.boro);
if ~isequal(size_stack, size_boro)
    errordlg('Las dimensiones de IMAGEN y DOSIS no coinciden.', ...
        'Error de dimensiones');
end

if sum(patient.D1c.boro(:))==0 
   errordlg('La matriz es cero' , ...
        'Error');
end 


% --- Executes on button press in btn_showMCNPerror.
function btn_showMCNPerror_Callback(hObject, eventdata, handles)
% hObject    handle to btn_showMCNPerror (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in btn_abs_dose.
global patient;

set(handles.txt_status,'String','Showing MCNP Boron Dose RelError...')
set(handles.rbtn_showMCNPerror,'Enable','off')
set(handles.rbtn_showMCNPerror,'Value',1)

% prompt = {'Error'};
% dlg_title = 'Error';
% num_lines = 1;
% defaultans = {'1'};
% answ = inputdlg(prompt,dlg_title,num_lines,defaultans);
% nE = str2double(answ{1});

list={'Boro','Thermal','Fast','Gamma'};
[nE,tf] = listdlg('PromptString','Error','ListString',list,'SelectionMode','single');

% patient.Dt.E(patient.stack == 0) = NaN;
[patient.xdmax,patient.ydmax,patient.zdmax] = ind2sub(size(patient.Dt.E(:,:,:,nE)),find(patient.Dt.E(:,:,:,nE)==min(min(min(patient.Dt.E(:,:,:,nE))))));


patient.px = patient.xdmax;
patient.py = patient.ydmax;
patient.nimg = patient.zdmax;
patient.nE=nE;

patient.show = 4;

if(isfield(patient,'himg'))
    patient = rmfield(patient,'himg');
end
if(isfield(patient,'dose_surf'))
    patient = rmfield(patient,'dose_surf');
end
if(isfield(patient,'hisoc'))
    patient = rmfield(patient,'hisoc');
end
if(isfield(patient,'hplotmax'))
    patient = rmfield(patient,'hplotmax');
end

ax=gca;
cla(ax)
ax.Toolbar.Visible = 'on';

plot_choosen_dose(hObject,1)

if(isfield(patient,'Dw'))
    set(handles.rbtn_showDw,'Enable','off')
    set(handles.rbtn_showDw,'Value',0)
end

if(isfield(patient,'Diso'))
    set(handles.rbtn_showDiso,'Enable','off')
    set(handles.rbtn_showDiso,'Value',0)
end

if(isfield(patient,'Df'))
    set(handles.rbtn_showDf,'Enable','off')
    set(handles.rbtn_showDf,'Value',0)
end

set(handles.txt_status,'String','Showing MCNP Boron Dose RelError...')
set(handles.btn_create_dvh,'Enable','off')
set(handles.btn_utcp,'Enable','off')
set(handles.btn_tir_point,'Enable','off')
set(handles.btn_find,'Enable','off')


% --- Executes on button press in rbtn_showMCNPerror.
function rbtn_showMCNPerror_Callback(hObject, eventdata, handles)
% hObject    handle to rbtn_showMCNPerror (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient;

set(handles.rbtn_showMCNPerror,'Enable','off')
set(handles.rbtn_showMCNPerror,'Value',1)

[patient.xdmax,patient.ydmax,patient.zdmax] = ind2sub(size(patient.Dt.E),find(patient.Dt.E==min(min(min(patient.Dt.E)))));

%patient.show = 4;
%plot_choosen_dose(hObject,1)

if(isfield(patient,'Dw'))
    set(handles.rbtn_showDw,'Enable','off')
    set(handles.rbtn_showDw,'Value',0)
end

if(isfield(patient,'Diso'))
    set(handles.rbtn_showDiso,'Enable','off')
    set(handles.rbtn_showDiso,'Value',0)
end

if(isfield(patient,'Df'))
    set(handles.rbtn_showDf,'Enable','off')
    set(handles.rbtn_showDf,'Value',0)
end

set(handles.txt_status,'String','Calculating Absorbed Dose...OK')
set(handles.btn_create_dvh,'Enable','off')
set(handles.btn_utcp,'Enable','off')
set(handles.btn_tir_point,'Enable','off')
set(handles.btn_find,'Enable','off')
% Hint: get(hObject,'Value') returns toggle state of rbtn_showMCNPerror

function btn_abs_dose_Callback(hObject, eventdata, handles)
% hObject    handle to btn_abs_dose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient;
set(handles.txt_status,'String','Calculating Absorbed Dose...')

patient.Dt.boro = patient.D1.boro*patient.boron1;


prompt = {'Tissue/Blood Ratio:'};
dlg_title = 'Absorbed Dose';
num_lines = 1;
defaultans = {'1'};
answ = inputdlg(prompt,dlg_title,num_lines,defaultans);
ratio = str2double(answ{1});

% Calculo Dosis Fisica.
patient.Df = ((ratio*patient.Dt.boro) + patient.Dt.fast + patient.Dt.thn + patient.Dt.g)*patient.tir;
% patient.Df(patient.stack == 0) = NaN;
% Calculo de las coordenadas del maximo de la dosis Fisica.
[patient.xdmax,patient.ydmax,patient.zdmax] = ind2sub(size(patient.Df),find(patient.Df==max(max(max(patient.Df)))));
patient.Dfmax = max(max(max(patient.Df)));
patient.px = patient.xdmax(1);
patient.py = patient.ydmax(1);
patient.nimg = patient.zdmax(1);
% Calculo NTCP.
% Calculo las proporciones de las componentes respecto al punto elegido.




patient.rb = ratio*patient.Dt.boro*patient.tir./patient.Df;
patient.rf = patient.Dt.fast*patient.tir./patient.Df;
patient.rt = patient.Dt.thn*patient.tir./patient.Df;
patient.rg = patient.Dt.g*patient.tir./patient.Df;

% Calculo NTCP
parametros_prescription = model_selector(4);
%[patient.Daf,patient.ntcp] = ntcp(patient.Df,patient.rb,patient.rf+patient.rt,patient.rg,parametros_ntcp);
patient.Deq = isodose(patient.Dt,patient.tir,ratio,parametros_prescription(1:14));



seleccion = get(handles.tbl_roi,'Data');
seleccion=seleccion(:,3);
%seleccion=cell2mat(seleccion);
seleccion = cellfun(@(x) double(x), seleccion);
if sum(seleccion)==0
    warndlg('You must select an ROI before proceeding.', 'Advertencia');
    return
end


patient.dvh_roi=seleccion;
rois = patient.roi_level(patient.dvh_roi == 1);
mask = nan(size(patient.stack));
rois_names = patient.roi_list(patient.dvh_roi == 1);


mask = nan(size(patient.stack));
for i=1:length(rois)
    mask(patient.stack==rois(i))=1;
end
Dmask = patient.Deq .* mask;
dvh_dose = Dmask(mask==1);
dv = patient.xstep * patient.ystep * patient.zstep;
roi_vol = dv*length(dvh_dose);

[patient.Def,patient.ntcp] = ntcp(dvh_dose, dv/roi_vol, parametros_prescription(15:18));

%% hay que pasarlo al DVH
%nNaN=sum(isnan(patient.Df(:)));
%if nNaN>0
%   warndlg('Se encontraron valores NaN en los datos', 'Advertencia');
%end


set(handles.rbtn_showDf,'Enable','off')
set(handles.rbtn_showDf,'Value',1)

patient.show = 1;
if(isfield(patient,'himg'))
    patient = rmfield(patient,'himg');
end
if(isfield(patient,'dose_surf'))
    patient = rmfield(patient,'dose_surf');
end
if(isfield(patient,'hisoc'))
    patient = rmfield(patient,'hisoc');
end
if(isfield(patient,'hplotmax'))
    patient = rmfield(patient,'hplotmax');
end
ax=gca;
cla(ax)
ax.Toolbar.Visible = 'on';

plot_choosen_dose(hObject,1)

if(isfield(patient,'Dt'))
    set(handles.rbtn_showMCNPerror,'Enable','off')
    set(handles.rbtn_showMCNPerror,'Value',0)
end

if(isfield(patient,'Dw'))
    set(handles.rbtn_showDw,'Enable','off')
    set(handles.rbtn_showDw,'Value',0)
end

if(isfield(patient,'Diso'))
    set(handles.rbtn_showDiso,'Enable','off')
    set(handles.rbtn_showDiso,'Value',0)
end

set(handles.txt_status,'String','Calculating Absorbed Dose...OK')
set(handles.btn_create_dvh,'Enable','on')
set(handles.btn_utcp,'Enable','on')
set(handles.btn_tir_point,'Enable','on')
set(handles.btn_find,'Enable','on')
% set(handles.bnt_tir_region,'Enable','on')

% --- Executes on button press in rbtn_showDf.
function rbtn_showDf_Callback(hObject, eventdata, handles)
global patient;
set(handles.rbtn_showDf,'Enable','off')
set(handles.rbtn_showDf,'Value',1)

%patient.show = 1;
%plot_choosen_dose(hObject,1)

if(isfield(patient,'Dt'))
    set(handles.rbtn_showMCNPerror,'Enable','off')
    set(handles.rbtn_showMCNPerror,'Value',0)
end

if(isfield(patient,'Dw'))
    set(handles.rbtn_showDw,'Enable','off')
    set(handles.rbtn_showDw,'Value',0)
end

if(isfield(patient,'Diso'))
    set(handles.rbtn_showDiso,'Enable','off')
    set(handles.rbtn_showDiso,'Value',0)
end
% hObject    handle to rbtn_showDf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtn_showDf

% --- Executes on button press in btn_w_dose.
function btn_w_dose_Callback(hObject, eventdata, handles)
% hObject    handle to btn_w_dose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient;
set(handles.txt_status,'String','Calculating RBEw Dose...')

patient.Dt.boro = patient.D1.boro*patient.boron1;


prompt = {'Boron:','Thermal Neutron:','Fast Neutron:','Gamma:','Tissue/Blood Ratio:'};
dlg_title = 'RBE factors';
num_lines = 1;
defaultans = {'3.8','3.2','3.2','1','3.5'};
% Hanna
%defaultans = {'3.8','3.2','3.2','1','1'}; % la informacion del ratio esta en cada voxel
answ = inputdlg(prompt,dlg_title,num_lines,defaultans);

%RBE Factors.
rbe.boro = str2double(answ{1});
rbe.thn = str2double(answ{2});
rbe.fastn = str2double(answ{3});
rbe.g = str2double(answ{4});
%RBEw Dose.
patient.Dw = ((str2double(answ{5})*patient.Dt.boro*rbe.boro) + patient.Dt.fast*rbe.fastn + patient.Dt.thn*rbe.thn + patient.Dt.g*rbe.g)*patient.tir;
% Mascara
Dose=patient.Dw;
Dose(patient.stack==0)=0;
patient.Dw=Dose;
% patient.Dw(patient.stack == 0)= NaN;
[patient.xdmax,patient.ydmax,patient.zdmax] = ind2sub(size(patient.Dw),find(patient.Dw==max(max(max(patient.Dw)))));
patient.Dwmax = max(max(max(patient.Dw)));
patient.nimg = patient.zdmax;
patient.px = patient.xdmax;
patient.py = patient.ydmax;

set(handles.rbtn_showDw,'Enable','off')
set(handles.rbtn_showDw,'Value',1)

patient.show = 2;
if(isfield(patient,'himg'))
    patient = rmfield(patient,'himg');
end
if(isfield(patient,'dose_surf'))
    patient = rmfield(patient,'dose_surf');
end
if(isfield(patient,'hisoc'))
    patient = rmfield(patient,'hisoc');
end
if(isfield(patient,'hplotmax'))
    patient = rmfield(patient,'hplotmax');
end

ax=gca;
cla(ax)
ax.Toolbar.Visible = 'on';

plot_choosen_dose(hObject,1)

if(isfield(patient,'Dt'))
    set(handles.rbtn_showMCNPerror,'Enable','off')
    set(handles.rbtn_showMCNPerror,'Value',0)
end

if(isfield(patient,'Df'))
    set(handles.rbtn_showDf,'Enable','off')
    set(handles.rbtn_showDf,'Value',0)
end

if(isfield(patient,'Diso'))
    set(handles.rbtn_showDiso,'Enable','off')
    set(handles.rbtn_showDiso,'Value',0)
end
set(handles.txt_status,'String','Calculating RBEw Dose...OK')
set(handles.btn_create_dvh,'Enable','on')
% set(handles.btn_utcp,'Enable','on')
set(handles.btn_tir_point,'Enable','on')
set(handles.btn_find,'Enable','on')
% set(handles.bnt_tir_region,'Enable','on')

% --- Executes on button press in rbtn_showDw.
function rbtn_showDw_Callback(hObject, eventdata, handles)
global patient;

%patient.show = 2;
%plot_choosen_dose(hObject,1)

if(isfield(patient,'Dt'))
    set(handles.rbtn_showMCNPerror,'Enable','off')
    set(handles.rbtn_showMCNPerror,'Value',0)
end

if(isfield(patient,'Df'))
    set(handles.rbtn_showDf,'Enable','off')
    set(handles.rbtn_showDf,'Value',0)
end

if(isfield(patient,'Diso'))
    set(handles.rbtn_showDiso,'Enable','off')
    set(handles.rbtn_showDiso,'Value',0)
end
% hObject    handle to rbtn_showDw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtn_showDw

% --- Executes on button press in btn_isoe_dose.
function btn_isoe_dose_Callback(hObject, eventdata, handles)
% hObject    handle to btn_isoe_dose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient;
set(handles.txt_status,'String','Calculating Iso Effective Dose...')

prompt = {'Tumor/blood B10 ratio:'};
dlgtitle = 'Tumor/blood B10 ratio';
dims = [1 35];
definput = {'3.5'};
bratio = str2double(inputdlg(prompt,dlgtitle,dims,definput));

patient.Dt.boro = patient.D1.boro*patient.boron1;

% Calculo Dosis Isoefectiva.
parametros = model_selector(1);
patient.Diso = isodose(patient.Dt,patient.tir,bratio,parametros);
% patient.Diso(patient.stack == 0) = NaN;

% Calculo de las coordenadas del maximo de la dosis iso.
[patient.xdmax,patient.ydmax,patient.zdmax] = ind2sub(size(patient.Diso),find(patient.Diso==max(max(max(patient.Diso)))));
patient.Disomax = max(max(max(patient.Diso)));
patient.px = patient.xdmax;
patient.py = patient.ydmax;
patient.nimg = patient.zdmax;

set(handles.rbtn_showDiso,'Enable','off')
set(handles.rbtn_showDiso,'Value',1)

patient.show = 3;
if(isfield(patient,'himg'))
    patient = rmfield(patient,'himg');
end
if(isfield(patient,'dose_surf'))
    patient = rmfield(patient,'dose_surf');
end
if(isfield(patient,'hisoc'))
    patient = rmfield(patient,'hisoc');
end
if(isfield(patient,'hplotmax'))
    patient = rmfield(patient,'hplotmax');
end

ax=gca;
cla(ax)
ax.Toolbar.Visible = 'on';

plot_choosen_dose(hObject,1)

if(isfield(patient,'Dt'))
    set(handles.rbtn_showMCNPerror,'Enable','off')
    set(handles.rbtn_showMCNPerror,'Value',0)
end

if(isfield(patient,'Dw'))
    set(handles.rbtn_showDw,'Enable','off')
    set(handles.rbtn_showDw,'Value',0)
end

if(isfield(patient,'Df'))
    set(handles.rbtn_showDf,'Enable','off')
    set(handles.rbtn_showDf,'Value',0)
end

set(handles.txt_status,'String','Calculating Iso Effective Dose...OK')
set(handles.btn_create_dvh,'Enable','on')
% set(handles.btn_utcp,'Enable','on')
set(handles.btn_tir_point,'Enable','on')
set(handles.btn_find,'Enable','on')
% set(handles.bnt_tir_region,'Enable','on')

% --- Executes on button press in rbtn_showDiso.
function rbtn_showDiso_Callback(hObject, eventdata, handles)
global patient;
%patient.show = 3;
%plot_choosen_dose(hObject,1)

if(isfield(patient,'Dt'))
    set(handles.rbtn_showMCNPerror,'Enable','off')
    set(handles.rbtn_showMCNPerror,'Value',0)
end

if(isfield(patient,'Dw'))
    set(handles.rbtn_showDw,'Enable','off')
    set(handles.rbtn_showDw,'Value',0)
end

if(isfield(patient,'Df'))
    set(handles.rbtn_showDf,'Enable','off')
    set(handles.rbtn_showDf,'Value',0)
end
% hObject    handle to rbtn_showDiso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtn_showDiso

% --- Executes when entered data in editable cell(s) in tbl_roi.
function tbl_roi_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to tbl_roi (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global patient;

switch eventdata.Indices(2)
    case 2
        patient.roi_list(eventdata.Indices(1)) = string(eventdata.EditData);
    case 3
        patient.dvh_roi(eventdata.Indices(1)) = double(eventdata.EditData);
    case 4
        patient.roi_tumor(eventdata.Indices(1)) = double(eventdata.EditData);
end

data = get(handles.tbl_roi,'data');
data(eventdata.Indices(1),eventdata.Indices(2)) = {eventdata.EditData};
set(handles.tbl_roi,'data',data);

% --- Executes on button press in btn_create_dvh.
function btn_create_dvh_Callback(hObject, eventdata, handles)
global patient;

switch patient.show
    case 1
        Dose = patient.Df;
        type = 'Absorbed Dose';
        unit = ' Gy';
    case 2
        Dose = patient.Dw;
        type = 'RBEw Dose';
        unit = ' Gy (RBE)';
    case 3
        Dose = patient.Diso;
        type = 'IsoEffective Dose';
        unit = ' Gy (isoE)';
    otherwise
        warndlg('Select a dose type first.', 'Warning');
        return;
end

data = get(handles.tbl_roi,'Data');
seleccion = cellfun(@double, data(:,3));
if sum(seleccion) == 0
    warndlg('Select at least one ROI.', 'Warning');
    return;
end
patient.dvh_roi = seleccion;

sel_idx = find(seleccion == 1);
n_sel = length(sel_idx);
voxel_vol = patient.xstep * patient.ystep * patient.zstep;

roi_masks = cell(n_sel, 1);
roi_names = cell(n_sel, 1);
roi_colors = zeros(n_sel, 3);
roi_tumor = false(n_sel, 1);

for i = 1:n_sel
    idx = sel_idx(i);
    level = patient.roi_level(idx);
    roi_masks{i} = patient.stack == level;
    roi_names{i} = char(patient.roi_list(idx));
    roi_tumor(i) = logical(double(data{idx, 4}));
    roi_colors(i,:) = patient.roi_color(idx, :);
end

R = f_dvh(Dose, roi_masks, roi_names, roi_colors, voxel_vol, roi_tumor);

parametros_ntcp = model_selector(3);
parametros_tcp = model_selector(2);

hFigure = figure('Name','DVH','NumberTitle','off','Color','w');
ax = axes('Parent',hFigure);
hold(ax,'on');

h_plots = gobjects(n_sel, 1);
x_max = 0;
for i = 1:n_sel
    if isempty(R(i).dvh_dose)
        continue;
    end
    h_plots(i) = plot(ax, R(i).dvh_dose, R(i).dvh_vol, ...
        'Color', R(i).color, 'LineWidth', 2);
    x_max = max(x_max, max(R(i).dvh_dose));

    dv = voxel_vol;
    roi_vol = R(i).vol_cc;

    if R(i).tumor
        [tcp_val, peud_val] = TCPc(R(i).dvh_dose, dv, [parametros_tcp{2:11}]);
        R(i).TCP = tcp_val;
        R(i).PEUD = peud_val;
    else
        [deffect, ntcp_val] = ntcp(R(i).dvh_dose, dv/roi_vol, [parametros_ntcp{2:5}]);
        R(i).Deffect = deffect;
        R(i).NTCP = ntcp_val;
    end
end

legend_str = cell(n_sel, 1);
metrics_data = cell(n_sel, 9);
for i = 1:n_sel
    if isempty(R(i).dvh_dose)
        legend_str{i} = R(i).name;
        metrics_data(i,:) = {R(i).name, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN};
        continue;
    end
    legend_str{i} = R(i).name;
    metrics_data(i,:) = {
        R(i).name, ...
        sprintf('%.2f', R(i).vol_cc), ...
        sprintf('%.2f', R(i).Dmax), ...
        sprintf('%.2f', R(i).Dmean), ...
        sprintf('%.2f', R(i).D98), ...
        sprintf('%.2f', R(i).D95), ...
        sprintf('%.2f', R(i).D50), ...
        sprintf('%.2f', R(i).D2), ...
        sprintf('%.2f', R(i).V20) };
end

x_max = max(x_max, 1);
xlim(ax, [0 x_max * 1.05]);
ylim(ax, [0 102]);
xlabel(ax, sprintf('%s [%s]', type, unit));
ylabel(ax, 'Volume [%]');
grid(ax, 'on');
box(ax, 'on');
legend(ax, h_plots(h_plots ~= 0), legend_str(h_plots ~= 0), 'Location', 'best');

metrics_pos = [0.15 0.02 0.7 0.12];
uitable('Parent', hFigure, 'Units', 'normalized', 'Position', metrics_pos, ...
    'Data', metrics_data, ...
    'ColumnName', {'ROI','Vol[cc]','Dmax','Dmean','D98','D95','D50','D2','V20'}, ...
    'ColumnWidth', {80,60,60,60,60,60,60,60,60}, ...
    'RowName', []);

hold(ax, 'off');

% Prescription calculation using the first ROI curve
valid_idx = find(h_plots ~= 0, 1);
if ~isempty(valid_idx)
    h = h_plots(valid_idx);
else
    h = gobjects(0);
end

msg = "Do you want to calculate the irradiation time?";
resp = questdlg(msg, 'Calcular tir', 'Yes', 'No', 'No');
if strcmp(resp, 'Yes') && isgraphics(h)
    tir = f_prescripcion(h);
    patient.tir = tir;
    set(handles.inp_tir, 'String', num2str(tir));
end

function inp_tir_Callback(hObject, eventdata, handles)
% hObject    handle to inp_tir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inp_tir as text
%        str2double(get(hObject,'String')) returns contents of inp_tir as a double
global patient;
patient.tir=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function inp_tir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inp_tir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in bnt_tir_region.
function bnt_tir_region_Callback(hObject, eventdata, handles)
% hObject    handle to bnt_tir_region (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in btn_tir_point.
function btn_tir_point_Callback(hObject, eventdata, handles)
% hObject    handle to btn_tir_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient;

switch patient.show
    case 1
        % Calculo las proporciones de las componentes respecto al punto elegido.
        Dp = patient.Df(patient.px,patient.py,patient.nimg);
        rb = patient.rb(patient.px,patient.py,patient.nimg);
        rf = patient.rf(patient.px,patient.py,patient.nimg);
        rt = patient.rt(patient.px,patient.py,patient.nimg);
        rg = patient.rg(patient.px,patient.py,patient.nimg);

        parametros_ntcp = model_selector(3);
        % Construct a questdlg with three options
        choice = questdlg('Dose prescription to...','Choose...','NTCP','Maximum','NTCP');
        % Handle response
        switch choice
            case 'NTCP'
                prompt = {'Enter NCTP [%] limit:'};
                dlg_title = 'NTCP limit';
                num_lines = 1;
                defaultans = {'5'};
                answ = inputdlg(prompt,dlg_title,num_lines,defaultans);
                lim = str2double(answ)/100;
                patient.tir = fzero(@(t) ntcp2(Dp*t,rb,rf+rt,rg,parametros_ntcp)-lim,60);
                set(handles.inp_tir,'String',num2str(patient.tir,'%0.2f'));
            case 'Maximum'
                prompt = {'Enter Maximum Absorbed Dose [Gy]:'};
                dlg_title = 'Maximum Absorbed Dose [Gy]:';
                num_lines = 1;
                defaultans = {'6'};
                answ = inputdlg(prompt,dlg_title,num_lines,defaultans);
                lim = str2double(answ);
                patient.tir = fzero(@(t) lim-(Dp*t),60);
                set(handles.inp_tir,'String',num2str(patient.tir,'%0.2f'));
        end
end

% --- Executes on button press in btn_load_2fr.
function btn_load_2fr_Callback(hObject, eventdata, handles)
% hObject    handle to btn_load_2fr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient;
set(handles.txt_status,'String','Loading Dose Matrix...')

% Handle response
if (~isfield(patient,'D2'))
    choice = questdlg('Select type of Dose Matrix:','Dose Matrix','Sera','Multicell','Sera');
    switch choice
        case 'Sera'
            [filename, pathname] = uigetfile('*.dat','Select Sera Dose Matrix');
            if isequal(filename,0)||isequal(pathname,0)
            else
                prompt = {'Neutrons:','Gamma'};
                dlg_title = 'Fraction 2 - Scaling Factors:';
                num_lines = 1;
                %defaultans = {'15'};
                answ = inputdlg(prompt,dlg_title,num_lines);

                sf_boron = str2double(answ{1});
                sf_g = str2double(answ{2});

                patient.D2c = leer_sera(fullfile(pathname,filename),sf_g,sf_boron);
                patient.D2c.E = zeros(size(patient.D2c.boro));
            end
        case 'Multicell'
            prompt = {'Thermal:','Fast:','Gamma:'};
            dlg_title = 'Fraction 2 - Scaling Factors:';
            num_lines = 1;
            defaultans = {'0.97','1.05','1.07'};
            answ = inputdlg(prompt,dlg_title,num_lines,defaultans);

            sthn = str2double(answ{1});
            sfast = str2double(answ{2});
            sg = str2double(answ{3});

            patient.D2c = leer_multicell(sthn,sfast,sg);
    end

    prompt = {'Blood B10 Concentration:'};
    dlg_title = 'Fraction 2';
    num_lines = 1;
    defaultans = {'15'};
    answ = inputdlg(prompt,dlg_title,num_lines,defaultans);

    patient.boron2 = str2double(answ{1});

    set(handles.inp_2ndb10,'Enable','on');
    set(handles.inp_2ndb10,'String',num2str(patient.boron2));

    xmp = patient.xstep/2;
    ymp = patient.ystep/2;
    zmp = patient.zstep/2;

    x=patient.xmin+xmp:patient.xstep:patient.xmax-xmp;
    y=patient.ymin+ymp:patient.ystep:patient.ymax-ymp;
    z=patient.zmin+zmp:patient.zstep:patient.zmax-zmp;

    [X,Y,Z] = ndgrid(x,y,z);
    [Xc,Yc,Zc]=ndgrid(patient.D2c.xc,patient.D2c.yc,patient.D2c.zc);

    % Interpolo componente a componente.
    set(handles.txt_status,'String','Interpolating...')

    patient.D2.boro = interpn(Xc,Yc,Zc,patient.D2c.boro,X,Y,Z,'linear');
    patient.D2.fast = interpn(Xc,Yc,Zc,patient.D2c.fast,X,Y,Z,'linear');
    patient.D2.thn = interpn(Xc,Yc,Zc,patient.D2c.thn,X,Y,Z,'linear');
    patient.D2.g = interpn(Xc,Yc,Zc,patient.D2c.g,X,Y,Z,'linear');
    patient.D2.E = interpn(Xc,Yc,Zc,patient.D2c.E,X,Y,Z,'linear');
end

prompt = {'Fraction 1 weight/time[min]:','Fraction 2 weight/time[min]:'};
dlg_title = 'Fraction weights';
num_lines = 1;
defaultans = {'0.5','0.5'};
answ = inputdlg(prompt,dlg_title,num_lines,defaultans);

fin1 = str2double(answ{1});
fin2 = str2double(answ{2});

patient.f1 = fin1/(fin1+fin2);
patient.f2 = fin2/(fin1+fin2);

txt = [num2str(patient.f1,'%0.2f') '*F1+' num2str(patient.f2,'%0.2f') '*F2'];
set(hObject,'String',txt);

patient.Dt.boro = (patient.f1*patient.D1.boro*patient.boron1)+(patient.f2*patient.D2.boro*patient.boron2);
patient.Dt.fast = (patient.f1*patient.D1.fast)+(patient.f2*patient.D2.fast);
patient.Dt.thn = (patient.f1*patient.D1.thn)+(patient.f2*patient.D2.thn);
patient.Dt.g = (patient.f1*patient.D1.g)+(patient.f2*patient.D2.g);
patient.Dt.E = patient.D2.E;

set(handles.inp_tir,'Enable','on')
set(handles.btn_showMCNPerror,'String','2nd BNCT Boron Relative Error')
set(handles.txt_status,'String','Loading Dose Matrix...OK')

function inp_xcor_Callback(hObject, eventdata, handles)
% hObject    handle to inp_xcor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inp_xcor as text
%        str2double(get(hObject,'String')) returns contents of inp_xcor as a double

% --- Executes during object creation, after setting all properties.
function inp_xcor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inp_xcor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLABintp
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function inp_ycor_Callback(hObject, eventdata, handles)
% hObject    handle to inp_ycor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of inp_ycor as text
%        str2double(get(hObject,'String')) returns contents of inp_ycor as a double

% --- Executes during object creation, after setting all properties.
function inp_ycor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inp_ycor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function inp_zcor_Callback(hObject, eventdata, handles)
% hObject    handle to inp_zcor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of inp_zcor as text
%        str2double(get(hObject,'String')) returns contents of inp_zcor as a double

% --- Executes during object creation, after setting all properties.
function inp_zcor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inp_zcor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in btn_go.
function btn_go_Callback(hObject, eventdata, handles)
% hObject    handle to btn_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient;
patient.px = str2double(get(handles.inp_xcor,'String'));
patient.py = str2double(get(handles.inp_ycor,'String'));
patient.nimg = str2double(get(handles.inp_zcor,'String'));
plot_choosen_dose(hObject,1)

function btn_find_Callback(hObject, eventdata, handles)
% --- Executes on button press in btn_find.
% hObject    handle to btn_find (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient;

switch patient.show
    case 1
        % Absorbed Dose.
        Dose = patient.Df;
        type = 'Absorbed Dose';
        unit = ' Gy';
        d = 'D=';
    case 2
        % RBEw Dose.
        Dose = patient.Dw;
        type = 'RBEw Dose';
        unit = ' Gy (RBE)';
        d = 'Dw=';
    case 3
        % IsoEffective Dose.
        Dose = patient.Diso;
        type = 'IsoEffective Dose';
        unit = ' Gy (isoE)';
        d = 'Diso=';
end

rois = patient.roi_level(patient.dvh_roi == 1);
if(rois>0)
    mask = nan(size(patient.stack));
    %     rois_names = patient.roi_list(patient.dvh_roi == 1);
    for i=1:length(rois)
        mask(patient.stack==rois(i))=1;
    end
else
    mask = ones(size(patient.stack));
end
Dose = Dose .* mask;

if(patient.show==1)
    prompt = {[type '[' unit ']'],'boron fraction:','neutron fraction:','gamma fraction:'};
    dlg_title = 'Find...';
    num_lines = 1;
    answ = inputdlg(prompt,dlg_title,num_lines);

    find_val = str2double(answ{1});
    find_fb = str2double(answ{2});
    find_fn = str2double(answ{3});
    find_fg = str2double(answ{4});


    txt = [d num2str(find_val,'%0.1f')];
    set(hObject,'String',txt);

    rb = patient.rb .* mask;
    rn = (patient.rf+patient.rt) .* mask;
    rg = patient.rg .* mask;
    %patient.find_matrix = abs(Dose-find_val);
    patient.find_matrix = sqrt(((Dose-find_val).^2)+((rb-find_fb).^2)+((rn-find_fn).^2)+((rg-find_fg).^2));
    patient.find_list = Dose(mask==1);
    patient.find_list = sort(reshape(patient.find_matrix,[1 numel(patient.find_matrix)]),'ascend');
    patient.find_ind = 1;
    [x,y,z] = ind2sub(size(patient.find_matrix),find(patient.find_matrix==patient.find_list(patient.find_ind )));

    set(handles.inp_xcor,'String',num2str(x));
    set(handles.inp_ycor,'String',num2str(y));
    set(handles.inp_zcor,'String',num2str(z));

    set(handles.btn_findl,'Enable','off');
    set(handles.btn_findr,'Enable','on');
elseif(patient.show~=1)
    [x,y,z] = ind2sub(size(Dose),find(Dose==max(max(max(Dose)))));
    set(handles.inp_xcor,'String',num2str(x));
    set(handles.inp_ycor,'String',num2str(y));
    set(handles.inp_zcor,'String',num2str(z));
end

% --- Executes on button press in btn_findl.
function btn_findl_Callback(hObject, eventdata, handles)
% hObject    handle to btn_findl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient;

if(patient.find_ind==length(patient.find_list))
    set(handles.btn_findr,'Enable','on');
end

patient.find_ind = patient.find_ind - 1;

[x,y,z] = ind2sub(size(patient.find_matrix),find(patient.find_matrix==patient.find_list(patient.find_ind )));

set(handles.inp_xcor,'String',num2str(x));
set(handles.inp_ycor,'String',num2str(y));
set(handles.inp_zcor,'String',num2str(z));

if(patient.find_ind==1)
    set(hObject,'Enable','off');
end

% --- Executes on button press in btn_findr.
function btn_findr_Callback(hObject, eventdata, handles)
% hObject    handle to btn_findr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient;

if(patient.find_ind==1)
    set(handles.btn_findl,'Enable','on');
end

patient.find_ind = patient.find_ind + 1;

[x,y,z] = ind2sub(size(patient.find_matrix),find(patient.find_matrix==patient.find_list(patient.find_ind )));

set(handles.inp_xcor,'String',num2str(x));
set(handles.inp_ycor,'String',num2str(y));
set(handles.inp_zcor,'String',num2str(z));

if(patient.find_ind==length(patient.find_list))
    set(hObject,'Enable','off');
end

% --- Executes on button press in btn_utcp.
function btn_utcp_Callback(hObject, eventdata, handles)
% hObject    handle to btn_utcp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient;

prompt = {'Total time [min]:','step [min]:'};
dlg_title = 'UTCP';
num_lines = 1;
defaultans = {'400','4'};
answ = inputdlg(prompt,dlg_title,num_lines,defaultans);

rois = patient.roi_level(patient.dvh_roi == 1);
mask = nan(size(patient.stack));

for i=1:length(rois)
    mask(patient.stack==rois(i))=1;
end

tiempo =0:str2num(answ{2}):str2num(answ{1});
tcp_t = zeros(1,length(tiempo));
ntcp_t = zeros(1,length(tiempo));

h=waitbar(0,'Calculando UTCP...');

% param = [2.468388470867692 0.000100508105922 3.090888558919538 0.000091579941911];
% TD50g = 39.8;
% mg = 0.17;

xm = patient.px;
ym = patient.py;
zm = patient.nimg;
dv = patient.xstep * patient.ystep * patient.zstep;

parametros_iso = model_selector(1);

prompt = {'Tumor/blood B10 ratio:'};
dlgtitle = 'Tumor/blood B10 ratio';
dims = [1 35];
definput = {'3.5'};
bratio = inputdlg(prompt,dlgtitle,dims,definput);

parametros_tcp = model_selector(2);
parametros_ntcp = model_selector(3);

for i=1:length(tiempo)

    Dmask2 = isodose(patient.Dt,tiempo(i),bratio,parametros_iso) .* mask;
    dvh_dose2 = Dmask2(mask==1);
    tcp_t(i) = TCPc(dvh_dose2,dv,parametros_tcp);
    %[~,ntcp_t(i)] = ntcp(patient.Df(xm,ym,zm)*(tiempo(i)*patient.f1),patient.rb(xm,ym,zm),patient.rf(xm,ym,zm)+patient.rt(xm,ym,zm),patient.rg(xm,ym,zm),parametros_ntcp);
    ntcp_t(i) = ntcp2(patient.Df(xm,ym,zm)*(tiempo(i)*patient.f1),patient.rb(xm,ym,zm),patient.rf(xm,ym,zm)+patient.rt(xm,ym,zm),patient.rg(xm,ym,zm),parametros_ntcp);
    waitbar(i/length(tiempo),h);
end

delete(h);
utcp_t=tcp_t.*(1-ntcp_t);
tmax = ind2sub(size(utcp_t),find(utcp_t==max(utcp_t)));
ob(1) = tiempo(tmax);
ob(2) = utcp_t(tmax);
ob(3) = tcp_t(tmax);
ob(4) = ntcp_t(tmax);

figure
plot(tiempo,tcp_t,'k')
hold on
plot(tiempo,ntcp_t,'k--')
plot(tiempo,utcp_t,'k','LineWidth',1.5)
plot([ob(1) ob(1)],[0 1],'k:')
plot(ob(1),ob(2),'k.')
plot(ob(1),ob(3),'k.')
plot(ob(1),ob(4),'k.')

legend('\it{TCP }','\it{NTCP }','\it{UTCP }');
xlabel('Treatment time [min]');
ylabel('Probability');
box on
grid on

function inp_1stb10_Callback(hObject, eventdata, handles)
% hObject    handle to inp_1stb10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inp_1stb10 as text
%        str2double(get(hObject,'String')) returns contents of inp_1stb10 as a double
global patient

set(handles.txt_status,'String','Updating boron content...')
patient.boron1 = str2double(get(hObject,'String'));

if (~isfield(patient,'D2'))
    patient.Dt.boro = patient.D1.boro*patient.boron1;
    patient.Dt.fast = patient.D1.fast;
    patient.Dt.thn = patient.D1.thn;
    patient.Dt.g = patient.D1.g;
else
    patient.Dt.boro = (patient.f1*patient.D1.boro*patient.boron1)+(patient.f2*patient.D2.boro*patient.boron2);
    patient.Dt.fast = (patient.f1*patient.D1.fast)+(patient.f2*patient.D2.fast);
    patient.Dt.thn = (patient.f1*patient.D1.thn)+(patient.f2*patient.D2.thn);
    patient.Dt.g = (patient.f1*patient.D1.g)+(patient.f2*patient.D2.g);
end

patient.show = 0;
plot_choosen_dose(hObject,1);

set(handles.txt_status,'String','Boron content updated...')


% --- Executes during object creation, after setting all properties.
function inp_1stb10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inp_1stb10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function inp_2ndb10_Callback(hObject, eventdata, handles)
% hObject    handle to inp_2ndb10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inp_2ndb10 as text
%        str2double(get(hObject,'String')) returns contents of inp_2ndb10 as a double
global patient

set(handles.txt_status,'String','Updating boron content...')
patient.boron2 = str2double(get(hObject,'String'));

patient.Dt.boro = (patient.f1*patient.D1.boro*patient.boron1)+(patient.f2*patient.D2.boro*patient.boron2);
patient.Dt.fast = (patient.f1*patient.D1.fast)+(patient.f2*patient.D2.fast);
patient.Dt.thn = (patient.f1*patient.D1.thn)+(patient.f2*patient.D2.thn);
patient.Dt.g = (patient.f1*patient.D1.g)+(patient.f2*patient.D2.g);

patient.show = 0;
plot_choosen_dose(hObject,1);

set(handles.txt_status,'String','Boron content updated...')


% --- Executes during object creation, after setting all properties.
function inp_2ndb10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inp_2ndb10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_mX.
function btn_mX_Callback(hObject, eventdata, handles)
% hObject    handle to btn_mX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient

set(handles.txt_status,'String','Mirroring X axis...')

patient.D1c.boro = patient.D1c.boro(end:-1:1,:,:);
patient.D1c.fast = patient.D1c.fast(end:-1:1,:,:);
patient.D1c.thn = patient.D1c.thn(end:-1:1,:,:);
patient.D1c.g = patient.D1c.g(end:-1:1,:,:);

if (isfield(patient,'D2'))
    patient.D2c.boro = patient.D2c.boro(end:-1:1,:,:);
    patient.D2c.fast = patient.D2c.fast(end:-1:1,:,:);
    patient.D2c.thn = patient.D2c.thn(end:-1:1,:,:);
    patient.D2c.g = patient.D2c.g(end:-1:1,:,:);
end

patient.D1.boro = patient.D1.boro(end:-1:1,:,:);
patient.D1.fast = patient.D1.fast(end:-1:1,:,:);
patient.D1.thn = patient.D1.thn(end:-1:1,:,:);
patient.D1.g = patient.D1.g(end:-1:1,:,:);


patient.Dt.boro = patient.Dt.boro(end:-1:1,:,:);
patient.Dt.fast = patient.Dt.fast(end:-1:1,:,:);
patient.Dt.thn = patient.Dt.thn(end:-1:1,:,:);
patient.Dt.g = patient.Dt.g(end:-1:1,:,:);

patient.show = 0;
plot_choosen_dose(hObject,1);
set(handles.txt_status,'String','X axis mirrored')
set(handles.btn_mX,'BackgroundColor', [0 1 0]);


% --- Executes on button press in btn_mZ.
function btn_mZ_Callback(hObject, eventdata, handles)
% hObject    handle to btn_mZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient

set(handles.txt_status,'String','Mirroring Z axis...')

patient.D1c.boro = patient.D1c.boro(:,:,end:-1:1);
patient.D1c.fast = patient.D1c.fast(:,:,end:-1:1);
patient.D1c.thn = patient.D1c.thn(:,:,end:-1:1);
patient.D1c.g = patient.D1c.g(:,:,end:-1:1);

if (isfield(patient,'D2'))
    patient.D2c.boro = patient.D2c.boro(:,:,end:-1:1);
    patient.D2c.fast = patient.D2c.fast(:,:,end:-1:1);
    patient.D2c.thn = patient.D2c.thn(:,:,end:-1:1);
    patient.D2c.g = patient.D2c.g(:,:,end:-1:1);
end

patient.D1.boro = patient.D1.boro(:,end:-1:1,:);
patient.D1.fast = patient.D1.fast(:,end:-1:1,:);
patient.D1.thn = patient.D1.thn(:,end:-1:1,:);
patient.D1.g = patient.D1.g(:,end:-1:1,:);

patient.Dt.boro = patient.Dt.boro(:,:,end:-1:1);
patient.Dt.fast = patient.Dt.fast(:,:,end:-1:1);
patient.Dt.thn = patient.Dt.thn(:,:,end:-1:1);
patient.Dt.g = patient.Dt.g(:,:,end:-1:1);

patient.show = 0;
plot_choosen_dose(hObject,1);
set(handles.txt_status,'String','Z axis mirrored')
set(handles.btn_mZ,'BackgroundColor', [0 1 0]);


% --- Executes on button press in btn_mY.
function btn_mY_Callback(hObject, eventdata, handles)
% hObject    handle to btn_mY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient

set(handles.txt_status,'String','Mirroring Y axis...')

patient.D1c.boro = patient.D1c.boro(:,end:-1:1,:);
patient.D1c.fast = patient.D1c.fast(:,end:-1:1,:);
patient.D1c.thn = patient.D1c.thn(:,end:-1:1,:);
patient.D1c.g = patient.D1c.g(:,end:-1:1,:);

if (isfield(patient,'D2'))
    patient.D2c.boro = patient.D2c.boro(:,end:-1:1,:);
    patient.D2c.fast = patient.D2c.fast(:,end:-1:1,:);
    patient.D2c.thn = patient.D2c.thn(:,end:-1:1,:);
    patient.D2c.g = patient.D2c.g(:,end:-1:1,:);
end

patient.D1.boro = patient.D1.boro(:,end:-1:1,:);
patient.D1.fast = patient.D1.fast(:,end:-1:1,:);
patient.D1.thn = patient.D1.thn(:,end:-1:1,:);
patient.D1.g = patient.D1.g(:,end:-1:1,:);


patient.Dt.boro = patient.Dt.boro(:,end:-1:1,:);
patient.Dt.fast = patient.Dt.fast(:,end:-1:1,:);
patient.Dt.thn = patient.Dt.thn(:,end:-1:1,:);
patient.Dt.g = patient.Dt.g(:,end:-1:1,:);

patient.show = 0;
plot_choosen_dose(hObject,1);
set(handles.txt_status,'String','Y axis mirrored')
set(handles.btn_mY,'BackgroundColor', [0 1 0]);

% --- Executes on button press in btn_sXY.
function btn_sXY_Callback(hObject, eventdata, handles)
% hObject    handle to btn_sXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient

set(handles.txt_status,'String','Switching X and Y axis...')

patient.D1c.boro = permute(patient.D1c.boro,[2 1 3]);
patient.D1c.fast = permute(patient.D1c.fast,[2 1 3]);
patient.D1c.thn = permute(patient.D1c.thn,[2 1 3]);
patient.D1c.g = permute(patient.D1c.g,[2 1 3]);

xmp = patient.xstep/2;
ymp = patient.ystep/2;
zmp = patient.zstep/2;

x=patient.xmin+xmp:patient.xstep:patient.xmax-xmp;
y=patient.ymin+ymp:patient.ystep:patient.ymax-ymp;
z=patient.zmin+zmp:patient.zstep:patient.zmax-zmp;

[X,Y,Z] = ndgrid(x,y,z);
[Xc,Yc,Zc]=ndgrid(patient.D1c.xc,patient.D1c.yc,patient.D1c.zc);

% Interpolo componente a componente.
set(handles.txt_status,'String','Interpolating 1st BNCT...')

patient.D1.boro = interpn(Xc,Yc,Zc,patient.D1c.boro,X,Y,Z,'linear');
patient.D1.fast = interpn(Xc,Yc,Zc,patient.D1c.fast,X,Y,Z,'linear');
patient.D1.thn = interpn(Xc,Yc,Zc,patient.D1c.thn,X,Y,Z,'linear');
patient.D1.g = interpn(Xc,Yc,Zc,patient.D1c.g,X,Y,Z,'linear');

if(isfield(patient,'D2'))
    patient.D2c.boro = permute(patient.D2c.boro,[2 1 3]);
    patient.D2c.fast = permute(patient.D2c.fast,[2 1 3]);
    patient.D2c.thn = permute(patient.D2c.thn,[2 1 3]);
    patient.D2c.g = permute(patient.D2c.g,[2 1 3]);

    set(handles.txt_status,'String','Interpolating 2nd BNCT...')

    patient.D2.boro = interpn(Xc,Yc,Zc,patient.D2c.boro,X,Y,Z,'linear');
    patient.D2.fast = interpn(Xc,Yc,Zc,patient.D2c.fast,X,Y,Z,'linear');
    patient.D2.thn = interpn(Xc,Yc,Zc,patient.D2c.thn,X,Y,Z,'linear');
    patient.D2.g = interpn(Xc,Yc,Zc,patient.D2c.g,X,Y,Z,'linear');

    patient.Dt.boro = (patient.f1*patient.D1.boro*patient.boron1)+(patient.f2*patient.D2.boro*patient.boron2);
    patient.Dt.fast = (patient.f1*patient.D1.fast)+(patient.f2*patient.D2.fast);
    patient.Dt.thn = (patient.f1*patient.D1.thn)+(patient.f2*patient.D2.thn);
    patient.Dt.g = (patient.f1*patient.D1.g)+(patient.f2*patient.D2.g);
else
    patient.Dt.boro = patient.D1.boro*patient.boron1;
    patient.Dt.fast = patient.D1.fast;
    patient.Dt.thn = patient.D1.thn;
    patient.Dt.g = patient.D1.g;
end
patient.tir = 1;
set(handles.inp_tir,'Enable','on')

patient.show = 0;
plot_choosen_dose(hObject,1);
set(handles.txt_status,'String','X and Y axis switched...')


% --- Executes on button press in btn_sZX.
function btn_sZX_Callback(hObject, eventdata, handles)
% hObject    handle to btn_sZX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient

set(handles.txt_status,'String','Switching X and Z axis...')

patient.D1c.boro = permute(patient.D1c.boro,[3 2 1]);
patient.D1c.fast = permute(patient.D1c.fast,[3 2 1]);
patient.D1c.thn = permute(patient.D1c.thn,[3 2 1]);
patient.D1c.g = permute(patient.D1c.g,[3 2 1]);

xmp = patient.xstep/2;
ymp = patient.ystep/2;
zmp = patient.zstep/2;

x=patient.xmin+xmp:patient.xstep:patient.xmax-xmp;
y=patient.ymin+ymp:patient.ystep:patient.ymax-ymp;
z=patient.zmin+zmp:patient.zstep:patient.zmax-zmp;

[X,Y,Z] = ndgrid(x,y,z);
[Xc,Yc,Zc]=ndgrid(patient.D1c.xc,patient.D1c.yc,patient.D1c.zc);

% Interpolo componente a componente.
set(handles.txt_status,'String','Interpolating 1st BNCT...')

patient.D1.boro = interpn(Xc,Yc,Zc,patient.D1c.boro,X,Y,Z,'linear');
patient.D1.fast = interpn(Xc,Yc,Zc,patient.D1c.fast,X,Y,Z,'linear');
patient.D1.thn = interpn(Xc,Yc,Zc,patient.D1c.thn,X,Y,Z,'linear');
patient.D1.g = interpn(Xc,Yc,Zc,patient.D1c.g,X,Y,Z,'linear');

if(isfield(patient,'D2'))
    patient.D2c.boro = permute(patient.D2c.boro,[3 2 1]);
    patient.D2c.fast = permute(patient.D2c.fast,[3 2 1]);
    patient.D2c.thn = permute(patient.D2c.thn,[3 2 1]);
    patient.D2c.g = permute(patient.D2c.g,[3 2 1]);

    set(handles.txt_status,'String','Interpolating 2nd BNCT...')

    patient.D2.boro = interpn(Xc,Yc,Zc,patient.D2c.boro,X,Y,Z,'linear');
    patient.D2.fast = interpn(Xc,Yc,Zc,patient.D2c.fast,X,Y,Z,'linear');
    patient.D2.thn = interpn(Xc,Yc,Zc,patient.D2c.thn,X,Y,Z,'linear');
    patient.D2.g = interpn(Xc,Yc,Zc,patient.D2c.g,X,Y,Z,'linear');

    patient.Dt.boro = (patient.f1*patient.D1.boro*patient.boron1)+(patient.f2*patient.D2.boro*patient.boron2);
    patient.Dt.fast = (patient.f1*patient.D1.fast)+(patient.f2*patient.D2.fast);
    patient.Dt.thn = (patient.f1*patient.D1.thn)+(patient.f2*patient.D2.thn);
    patient.Dt.g = (patient.f1*patient.D1.g)+(patient.f2*patient.D2.g);
else
    patient.Dt.boro = patient.D1.boro*patient.boron1;
    patient.Dt.fast = patient.D1.fast;
    patient.Dt.thn = patient.D1.thn;
    patient.Dt.g = patient.D1.g;
end
patient.tir = 1;
set(handles.inp_tir,'Enable','on')

patient.show = 0;
plot_choosen_dose(hObject,1);
set(handles.txt_status,'String','X and Z axis switched...')


% --- Executes on button press in btn_sYZ.
function btn_sYZ_Callback(hObject, eventdata, handles)
% hObject    handle to btn_sYZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient

set(handles.txt_status,'String','Switching Y and Z axis...')

patient.D1c.boro = permute(patient.D1c.boro,[1 3 2]);
patient.D1c.fast = permute(patient.D1c.fast,[1 3 2]);
patient.D1c.thn = permute(patient.D1c.thn,[1 3 2]);
patient.D1c.g = permute(patient.D1c.g,[1 3 2]);

xmp = patient.xstep/2;
ymp = patient.ystep/2;
zmp = patient.zstep/2;

x=patient.xmin+xmp:patient.xstep:patient.xmax-xmp;
y=patient.ymin+ymp:patient.ystep:patient.ymax-ymp;
z=patient.zmin+zmp:patient.zstep:patient.zmax-zmp;

[X,Y,Z] = ndgrid(x,y,z);
[Xc,Yc,Zc]=ndgrid(patient.D1c.xc,patient.D1c.yc,patient.D1c.zc);

% Interpolo componente a componente.
set(handles.txt_status,'String','Interpolating 1st BNCT...')

patient.D1.boro = interpn(Xc,Yc,Zc,patient.D1c.boro,X,Y,Z,'linear');
patient.D1.fast = interpn(Xc,Yc,Zc,patient.D1c.fast,X,Y,Z,'linear');
patient.D1.thn = interpn(Xc,Yc,Zc,patient.D1c.thn,X,Y,Z,'linear');
patient.D1.g = interpn(Xc,Yc,Zc,patient.D1c.g,X,Y,Z,'linear');

if(isfield(patient,'D2'))
    patient.D2c.boro = permute(patient.D2c.boro,[1 3 2]);
    patient.D2c.fast = permute(patient.D2c.fast,[1 3 2]);
    patient.D2c.thn = permute(patient.D2c.thn,[1 3 2]);
    patient.D2c.g = permute(patient.D2c.g,[1 3 2]);

    set(handles.txt_status,'String','Interpolating 2nd BNCT...')

    patient.D2.boro = interpn(Xc,Yc,Zc,patient.D2c.boro,X,Y,Z,'linear');
    patient.D2.fast = interpn(Xc,Yc,Zc,patient.D2c.fast,X,Y,Z,'linear');
    patient.D2.thn = interpn(Xc,Yc,Zc,patient.D2c.thn,X,Y,Z,'linear');
    patient.D2.g = interpn(Xc,Yc,Zc,patient.D2c.g,X,Y,Z,'linear');

    patient.Dt.boro = (patient.f1*patient.D1.boro*patient.boron1)+(patient.f2*patient.D2.boro*patient.boron2);
    patient.Dt.fast = (patient.f1*patient.D1.fast)+(patient.f2*patient.D2.fast);
    patient.Dt.thn = (patient.f1*patient.D1.thn)+(patient.f2*patient.D2.thn);
    patient.Dt.g = (patient.f1*patient.D1.g)+(patient.f2*patient.D2.g);
else
    patient.Dt.boro = patient.D1.boro*patient.boron1;
    patient.Dt.fast = patient.D1.fast;
    patient.Dt.thn = patient.D1.thn;
    patient.Dt.g = patient.D1.g;
end
patient.tir = 1;
set(handles.inp_tir,'Enable','on')

patient.show = 0;
plot_choosen_dose(hObject,1);
set(handles.txt_status,'String','Y and Z axis switched...')


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global patient

if(exist('patient','var'))
    clear patient
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function img_CreateFcn(hObject, eventdata, handles)
% hObject    handle to img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate img
