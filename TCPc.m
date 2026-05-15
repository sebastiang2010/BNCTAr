%% TCP model for inhomogeneous dose distribution %%
% Mathematical model reported in G F Perotti Bernardini et al 2022 Radiat. Res. 198 

% Function to calculate tumor control probability (TCP) 
% & probability-equivalent uniform dose (PEUD)
% input: 
% % Deq: Tumor dose vector in photon equivalent units 
% % dv: voxel volume (cm^3) 
% % parameters: model parameters included in model_list.txt

function [TCP, PEUD]=TCPc(Deq,dv,parameters)

%% ***** Block of parameters ***** %% 
% alpha and beta for reference photons 
Ar = parameters(1);  
Br = parameters(2); 

% TCP volume parameters 
c1 = parameters(3);
c2 = parameters(4);

% Repair times & proportions for time factor
tof= parameters(5); % fast repair time
tos= parameters(6); % slow repair time 
% Low LET proportions
pf_g=parameters(7); 
ps_g=parameters(8);
% High LET proportions
pf_bnct=parameters(9);
ps_bnct=parameters(10);


%% ***** Lea-Catcheside Time factor ***** %%

%%%% Gr: photon reference time factor %%%%
% time (min): representative time of photon irradiation
time=30;

xf=tof./time;
xs=tos./time;                   
Gf=2*xf.*(1-xf.*(1-exp(-1./xf)));
Gs=2*xs.*(1-xs.*(1-exp(-1./xs)));

% contribution per component
f_g=1; 
f_bnct=0;

Gr=Gs-(f_g*pf_g+f_bnct*pf_bnct)*(Gs-Gf);


%% ***** TCP calculation ***** %%


Deq=Deq(~isnan(Deq));

H = exp(-Ar.*Deq-Gr.*Br.*(Deq.^2)).^(1/c2); 

TCP=exp(-c1*(sum(H.*dv)).^c2);

%% ***** PEUD calculation ***** %%
V=dv*length(Deq); % total tumor volume

C= log(-log(TCP)/(c1*V^c2));
PEUD = (-Ar + (Ar^2 - 4*Gr*Br*C)^0.5)/(2*Gr*Br);



