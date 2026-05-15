%% Photon isoeffective dose calculation %%
% Mathematical dose model reported in Gonzalez - Santa Cruz 2012 Radiat. Res. 178

% Function to calculate photon isoeffective doses (DisoE)
% input:
% % D: physical dose matrices (Gy) (struct) 
% % % D.boro: physical boron dose multiplied by blood boron or tissue boron (bratio=1) concentration ([concentration] = ppm) 
% % % D.thn: physical thermal neutron dose
% % % D.fast: physical fast neutron dose
% % % D.g: physical gamma dose 
% % tir: irradiation time (minutes)
% % bratio: tissue/blood boron concentration ratio 
% % parameters: model parameters included in model_list.txt

function DisoE = isodose(D,tir,bratio,parameters)

%% ***** Block of parameters ***** %% 
% alpha and beta for reference photons 
alpha_r = parameters(1);
beta_r = parameters(2);
alphabeta_r = alpha_r/beta_r;

% alpha and beta for BNCT components 
% (boron, thermal neutron: thn, fast neutron: nfast) 
alpha_boro = parameters(3);
beta_boro = parameters(4);
alpha_thn = parameters(5);
beta_thn = parameters(6);
alpha_fast = parameters(7);
beta_fast = parameters(8);

alpha = [alpha_boro,alpha_thn,alpha_fast,alpha_r,alpha_r];
beta = [beta_boro,beta_thn,beta_fast,beta_r,beta_r];

% SLD repair times & proportions for time factor
tof = parameters(9);  % fast repair time
tos = parameters(10); % slow repair time
% Low LET proportions
pf_g = parameters(11);
ps_g = parameters(12);
% High LET proportions
pf_bnct = parameters(13);
ps_bnct = parameters(14);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Gij: BNCT time factor %%%%

% 3= gamma (low LET)
% 2= neutron (thn + fast) (high LET)
% 1= boron (high LET)

xf=tof./tir;
xs=tos./tir;
Gf=2*xf.*(1-xf.*(1-exp(-1./xf)));
Gs=2*xs.*(1-xs.*(1-exp(-1./xs)));

% total physical dose
Dt = (D.boro*bratio)+D.thn+D.fast+D.g;
% dose contributions per bnct component
fb=(D.boro*bratio)./Dt; 
fn=(D.thn+D.fast)./Dt;
fg=D.g./Dt;

% relative proportions between components
a_12=fb./(fb+fn);
a_21=fn./(fb+fn);
a_13=fb./(fb+fg);
a_31=fg./(fb+fg);
a_23=fn./(fn+fg);
a_32=fg./(fn+fg);

% Gii
G_11=pf_bnct.*Gf+ps_bnct*Gs;
G_22=G_11;
G_33=pf_g.*Gf+ps_g*Gs;

% Gij
G_12=Gs-(a_12.*pf_bnct+a_21.*pf_bnct).*(Gs-Gf);
G_23=Gs-(a_23.*pf_bnct+a_32.*pf_g).*(Gs-Gf);
G_31=Gs-(a_31.*pf_g+a_13.*pf_bnct).*(Gs-Gf);

%% ***** Isoeffective dose calculation ***** %%

% linear term
Dosis_Tlineal = alpha(:,1)*(D.boro*bratio)+...
                alpha(:,2)*D.thn+...
                alpha(:,3)*D.fast+...
                alpha(:,4)*D.g;
Dosis_Tlineal=Dosis_Tlineal.*tir;

% quadratic term
Dosis_Tcuad=G_11.*beta(:,1).*((D.boro*bratio)).^2+...
    G_22.*beta(:,2).*(D.thn).^2+...
    G_22*beta(:,3).*(D.fast).^2+...
    G_33.*beta(:,4).*(D.g).^2+...
    2*G_12.*(beta(:,1)*beta(:,2)).^0.5.*((D.boro*bratio)).*(D.thn)+...
    2*G_12.*(beta(:,1)*beta(:,3)).^0.5.*((D.boro*bratio)).*(D.fast)+...
    2*G_31.*(beta(:,1)*beta(:,4)).^0.5.*((D.boro*bratio)).*(D.g)+...
    2*G_11.*(beta(:,2)*beta(:,3)).^0.5.*(D.thn).*(D.fast)+...
    2*G_23.*(beta(:,2)*beta(:,4)).^0.5.*(D.thn).*(D.g)+...
    2*G_23.*(beta(:,3)*beta(:,4)).^0.5.*(D.fast).*(D.g);
Dosis_Tcuad=Dosis_Tcuad.*tir.^2;

DisoE =0.5*(alphabeta_r./Gr).*((1+(4*Gr/alpha_r/alphabeta_r).*(Dosis_Tlineal+Dosis_Tcuad)).^0.5-1);

% nNaN=sum(isnan(DisoE(:))); 
% 
% if nNaN>0
%    warndlg('Se encontraron valores NaN en los datos', 'Advertencia');
% end 
end     
