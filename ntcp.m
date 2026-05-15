%% NTCP model %%
% Mathematical model reported in:
% S J González et al 2017 Phys. Med. Biol. 62 7938 (maximum dose) &
% A M Dattoli Viegas et al 2024 Physica Medica Phys. Med. 128 104840 (EUD)

% Function to calculate NTCP value (y) & prescription dose (D_presc)
% input: 
% % Deq: Normal tissue dose vector in photon equivalent units 
% % vi: volume fraction of voxel i with respect to total normal tissue volume (= voxel volume / total region volume)
% % parameters: model parameters included in model_list.txt

function [ D_presc , y ] = ntcp(Deq , vi, parameters)

%% ***** Block of parameters ***** %% 
% NTCP parameters 
TD50 = parameters(1);
m = parameters(2);

% alpha/beta for reference photons 
ab = parameters(3);

% EUD parameter (if the effect is asosiated with maximum dose a=0)
a = parameters(4);

%% ***** D_presc & NTCP calculation ***** %%

% NTCP model selection
if a == 0
% Mucosa NTCP model, prescription dose to the maximum
    D_presc = max(Deq(:));
    EQD2 = D_presc.*(ab+D_presc)./(ab+2);
    s=(EQD2-TD50)./(m*TD50);
else
% Brain NTCP model, prescription dose to the EUD
    D_presc = (sum(Deq.^a.*vi)).^(1/a);
    EUD = D_presc;
    s = (EUD - TD50)/(m*TD50)/2^0.5;
end

y=1./2*(1+erf(s));
end

