function [ y ] = ntcp2(Dtot,fb, fn, fg, parametros)
% UNTITLED3 Summary of this function goes here
% Detailed explanation goes here

tof = parametros(1);
tos = parametros(2);
pf_g = parametros(3);
ps_g = parametros(4);
pf_bnct = parametros(5);
ps_bnct = parametros(6);
param(1) = parametros(7);
param(2) = parametros(8);
param(3) = parametros(9);
param(4) = parametros(10);
TD50g = parametros(11);
mg = parametros(12);
% a y aj son el alpha R y el alpha/beta R del paper de Strigari.
a=0.35;
ab=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       CÁLCULO DEL Gr           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Time: tiempo referencia irrad con fotones
time=30;
xf=tof./time;
xs=tos./time;
                    
Gf=2*xf.*(1-xf.*(1-exp(-1./xf)));
Gs=2*xs.*(1-xs.*(1-exp(-1./xs)));
% Contribuciones de cada componente
p_gr=1;
p_bnctr=0;

% i=gamma bajo let
% j=bnct alto let
Gr=Gs-(p_gr*pf_g+p_bnctr*pf_bnct)*(Gs-Gf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       CÁLCULO DEL Gij           %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xf=tof./time;
xs=tos./time;
                    
Gf=2*xf.*(1-xf.*(1-exp(-1./xf)));
Gs=2*xs.*(1-xs.*(1-exp(-1./xs)));

% Contribuciones de cada componente
%fb, fn, fg
p_bnctB=fb;
p_bnctN=fn;
p_g=fg;

% 3=gamma bajo let
% 2=neutorn bnct alto let
% 1= boro alto let

% Gii
G_11=pf_bnct.*Gf+ps_bnct*Gs;
G_22=G_11;
G_33=pf_g.*Gf+ps_g*Gs;

G_12=Gs-(a_12.*pf_bnct+a_21.*pf_bnct).*(Gs-Gf);
G_23=Gs-(a_23.*pf_bnct+a_32.*pf_g).*(Gs-Gf);
G_31=Gs-(a_31.*pf_g+a_13.*pf_bnct).*(Gs-Gf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 
% to=46/0.693;  %3.2*60;
% time=40;
% p=to./time;
%                     
% G=2*p.*(1-p.*(1-exp(-1./p))); 
% Gr=0.95;

% param=[alpha_n beta_n^0.5 alpha_b beta_b^0.5]

% A=(fb*param(3)+fn*param(1)+fg*a);
% B=(fb*param(4)+fn*param(2)+fg*(a/10).^0.5).^2;
% B(B==0)=1e-10;


% param=[alpha_n 0 alpha_b 0] 
% A=(fb*param(2)+fn*param(1)+fg*a);
% B=(fg*(a/10));

% param=[alpha_n alpha_g] y beta_g=alpha_g/10;
% A=(fb*param(3)+fn*param(1)+fg*param(2));
% B=(fg*(param(2)/10));

% EQD2=(Dtot.*((A./B)+G*Dtot))./((A./B)+2*G);
% EQD2=(Dtot.*(A+G*B.*Dtot))./(A+2*G*B);

% Dsf=0.5*(ab/Gr)*((1+(4*Gr/ab/a).*((fb*param(3)+fn*param(1)+fg*a).*Dtot+G.*((fb*param(4).^2+fn*param(2).^2+fg*(a/10).^0.5).^2).*Dtot.^2)).^0.5-1);

A=(G_11.*fb.^2.*param(4).^4+G_22.*fn.^2.*param(2).^4+G_33.*fg.^2.*(a/10)+...
    2*G_12.*fb.*fn.*(param(4).*param(2)).^2+2*G_23.*fn.*fg.*(param(2).^4.*(a/10)).^0.5+...
    2*G_31.*fg.*fb.*((a/10).*param(4).^4).^0.5);

Dsf=0.5*(ab/Gr)*((1+(4*Gr/ab/a).*((fb*param(3)+fn*param(1)+fg*a).*Dtot+A.*Dtot.^2)).^0.5-1);



 Dloca=Dsf.*(ab+Dsf)./(ab+2);
s=(Dloca-TD50g)./(mg*TD50g);

% t=sum(isnan(s))
 
    
% s=(EQD2-TD50g)./(mg*TD50g);
y=1./2*(1+erf(s));
end

