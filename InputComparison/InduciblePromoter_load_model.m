%This script specifies the model MIP,r, having:
%   - structure M3D and 
%   - parameter estimates set to the values identified when fitting the structure on the pseudo-data from MPLac,r

model.input_model_type='charmodelC';                                        % Model introduction: 'charmodelC'|'c_model'|'charmodelM'|'matlabmodel'|'sbmlmodel'|                        'blackboxmodel'|'blackboxcost                             
model.n_st=3;                                                               % Number of states      
model.n_par=8;                                                              % Number of model parameters 
model.n_stimulus=1;                                                         % Number of inputs, stimuli or control variables   
model.st_names=char('Cit_mrna','Cit_foldedP','Cit_fluo');                   % Names of the states                                        
model.par_names=char('alpha1','Vm1','h1','Km1','d1','alpha2','d2','Kf');    % Names of the parameters                     
model.stimulus_names=char('IPTG');                                          % Names of the stimuli, inputs or controls                      
model.eqns=...                                                              % Equations describing system dynamics. Time derivatives are regarded 'd'st_name''
               char('dCit_mrna=alpha1+Vm1*(IPTG^h1/(Km1^h1+IPTG^h1))-d1*Cit_mrna',...
                    'dCit_foldedP=alpha2*Cit_mrna-(d2+Kf)*Cit_foldedP',...
                    'dCit_fluo=Kf*Cit_foldedP-d2*Cit_fluo');           
         
model.par = [0.0164186333380725 0.291556643109224 1.71763487775568 5.14394334860864 0.229999999999978 6.63776658557266 0.00575139649497780 0.0216999999961899]; % The parameter vector is set to the best estimates for MIP,r
