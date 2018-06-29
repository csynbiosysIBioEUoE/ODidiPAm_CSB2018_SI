function [ res ] = PLac_Compute_SteadyState( theta, IPTG )
%PLac_SteadyState Calculates the steady state for theta and IPTGe
%   Computes the steady state of the MPLac,r model for the given values of
%   theta and the input IPTG.

kLacI = theta(1);
k2 = theta(2);
kd = theta(3);
km2 = theta(4);
k1 = theta(5);
km1 = theta(6);
kLac12 = theta(7);
kTP1 = theta(8);
kcat = theta(9);
Km = theta(10);
kout = theta(11);
kC = theta(12);
lk = theta(13);
sc_molec = theta(14);

%%

K2 = km2/k2; 
K1 = km1/k1;

Lac12 = kLac12/(kTP1+kd); 
Lac12m = kTP1*Lac12/kd; 
IPTGi = kcat/(kout*kd)*Lac12m*IPTG/(IPTG+Km);
Ltot = kLacI/kd; 

L0 = Ltot/(1+(IPTGi/K2))^2; 
L1 = 2*L0*IPTGi/K2; 
L2 = (L0*(IPTGi)^2)/(K2^2); 

G20 = 1/((1+L0/K1)^2);
G21 = 2*L0*G20/K1;
G22 = ((L0/K1)^2)*G20;


Cit_molec = (kC*G20+lk*kC*(G21+G22))/kd; 
Cit_AU= sc_molec*Cit_molec; 

res = [L0 L1 L2 Lac12 Lac12m G20 G21 G22 IPTGi Cit_molec Cit_AU];
end

