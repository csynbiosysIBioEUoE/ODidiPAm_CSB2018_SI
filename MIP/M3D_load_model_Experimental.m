% As the experimental data are acquired in A.U., M3D is modified by adding
% an algebraic equation for the conversion from molecules number to A.U.
% This means that the parameter vector has an additional entry, represented
% by the scaling factor sc_molec

model.input_model_type='charmodelC';                                                % Model introduction: 'charmodelC'|'c_model'|'charmodelM'|'matlabmodel'|'sbmlmodel'|'blackboxmodel'|'blackboxcost                             
model.n_st=4;                                                                       % Number of states      
model.n_par=9;                                                                      % Number of model parameters 
model.n_stimulus=1;                                                                 % Number of inputs, stimuli or control variables   
model.st_names=char('Cit_mrna','Cit_foldedP','Cit_fluo','Cit_AU');                  % Names of the states                                        
model.par_names=char('alpha1','Vm1','h1','Km1','d1',...
                            'alpha2','d2','Kf','sc_molec');                         % Names of the parameters                     
model.stimulus_names=char('IPTG');                                                  % Names of the stimuli, inputs or controls                      
model.eqns=...                                                                      % Equations describing system dynamics. Time derivatives are regarded 'd'st_name''
               char('dCit_mrna=alpha1+Vm1*(IPTG^h1/(Km1^h1+IPTG^h1))-d1*Cit_mrna',...
                    'dCit_foldedP=alpha2*Cit_mrna-(d2+Kf)*Cit_foldedP',...
                    'dCit_fluo=Kf*Cit_foldedP-d2*Cit_fluo',...
                    'dCit_AU = sc_molec*dCit_fluo');           
                
%          alpha1	Vm1     h1  Km1     d1     alpha2  d2     Kf   sc_molec            
model.par =[0.000377125304442752, 0.00738924359598526, 1.53333782244337, 5.01927275636639, 0.00118831480244382, 0.0461264539194078, 0.000475563708997018, 0.000301803966012407, 68.8669567134881]; % results of the fitting of the reduced model performed using CrossValidation equal on Stelling Data
