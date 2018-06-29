function [ res ] = InduciblePromoter_steady_state( theta, IPTG )
%InduciblePromoter_steady_state Calculates the steady state for the theta and
%IPTG passed as inputs to the function

a1 = theta(1);
Vm1 = theta(2);
h1 = theta(3);
Km1 = theta(4);
d1 = theta(5);
a2 = theta(6);
d2 = theta(7);
Kf = theta(8);

% Algebraic equations are obtained by setting the time derivative of the
% state variables to 0.

Cit_mrna = (a1 + Vm1*(IPTG^h1/(Km1^h1+IPTG^h1)))/d1;

Cit_foldedP = (a2*Cit_mrna)/(Kf+d2);

Cit_fluo = (Kf*Cit_foldedP)/d2;


res = [Cit_mrna Cit_foldedP Cit_fluo];

end

