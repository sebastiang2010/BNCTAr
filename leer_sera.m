function Dsera = leer_sera(filename,gamma_cont,sf)
%% LEO DOSIMETRIA SERA.
% filename = '/Users/lucasprovenzano/Documents/Doctorado/Casos Finlandia/Head and neck/24_HN_11033M/1fr/r7_11a6_11_4753after/24_HN_11033M_r7_11a6_11_4753afterplan.dat';
% gamma_cont = 2.177;

%% CAMBIO A FORMATO DE COMPONENTES.
delimiter = ' ';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
data = [dataArray{1:end-1}];
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
patientdata = data(:,1:9);

% Data format of patientdata variable: #x y z gam hyd tot b10 gp n14 
patientdata_totphot = patientdata(:,4)+patientdata(:,8);

sca_b10 = sf;
sca_n14 = sf;
sca_hyd = sf;

% sca_b10 = 1.8609e13;
% sca_n14 = 1.8609e13;
% sca_hyd = 1.8609e13;

% melanoma
% sca_b10=1.35954E+13;
% sca_n14=1.35954E+13;
% sca_hyd=1.35954E+13;
 
% sca_b10=1.29480E+13;
% sca_n14=1.29480E+13;
% sca_hyd=1.29480E+13;

xdmax=max(patientdata(:,1));
xdmin=min(patientdata(:,1));

ydmax=max(patientdata(:,2));
ydmin=min(patientdata(:,2));

zdmax=max(patientdata(:,3));
zdmin=min(patientdata(:,3));

Dsera.xc=xdmin:xdmax;
Dsera.yc=ydmin:ydmax;
Dsera.zc=zdmin:zdmax;

ref_phot = patientdata_totphot(patientdata(:,6)==max(patientdata(:,6)));
sca_phot = gamma_cont/ref_phot;

Dsera.boro = permute(reshape(patientdata(:,7),30,30,30),[3 2 1])*sca_b10*1e-2;
Dsera.thn = permute(reshape(patientdata(:,9),30,30,30),[3 2 1])*sca_n14*1e-2;
Dsera.fast = permute(reshape(patientdata(:,5),30,30,30),[3 2 1])*sca_hyd*1e-2;
Dsera.g = permute(reshape(patientdata_totphot,30,30,30),[3 2 1])*sca_phot*1e-2;
