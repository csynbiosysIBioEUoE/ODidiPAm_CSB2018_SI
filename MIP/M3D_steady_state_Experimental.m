function [ res ] = M3D_steady_state_Experimental( theta, IPTG )
%   Computes the steady state of the M3D model structure for the given values of
%   theta and the input IPTG.

a1 = theta(1);
Vm1 = theta(2);
h1 = theta(3);
Km1 = theta(4);
d1 = theta(5);
a2 = theta(6);
d2 = theta(7);
Kf = theta(8);
sc_molec = theta(9);
                        
Cit_mrna = (a1 + Vm1*(IPTG^h1/(Km1^h1+IPTG^h1)))/d1;

Cit_foldedP = (a2*Cit_mrna)/(Kf+d2);

Cit_fluo = (Kf*Cit_foldedP)/d2;

Cit_AU= sc_molec*Cit_fluo; 

res = [Cit_mrna Cit_foldedP Cit_fluo Cit_AU];

end

