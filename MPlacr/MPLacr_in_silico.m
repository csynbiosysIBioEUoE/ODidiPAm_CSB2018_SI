%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to simulate the MPLacr model in response to a step in IPTG
% concentration specified in Run_MPLacr_in_silico_experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global epccOutputResultFileNameBase;
global IPTGe; 

resultFileName = [epccOutputResultFileNameBase,'.dat']
clear model;
clear exps;
clear sim;

results_folder = strcat('PLac_rep',datestr(now,'yyyy-mm-dd-HHMMSS'));
short_name     = 'PLac';

% Read the model into the model variable
PLac_load_model;

% Compile the model
clear inputs;
inputs.model = model;
inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = 'initial_setup';

AMIGO_Prep(inputs);
      
% Fixed parts of the experiment
duration = 24*60*60;     % Duration of the experiment in second

clear newExps;
newExps.n_exp = 1;
newExps.n_obs{1}=1;                                     % Number of observables per experiment                         
newExps.obs_names{1}=char('Cit_AU');                 % Name of the observables 
newExps.obs{1}=char('Fluorescence=Cit_AU');          % Observation function
newExps.exp_y0{1}=PLac_Compute_SteadyState(model.par, 0);    
    
newExps.t_f{1}=duration;                % Experiment duration
    
newExps.u_interp{1}='sustained';
newExps.u{1}=IPTGe;
newExps.t_con{1}=[0 duration]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mock an experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
clear inputs;
inputs.exps = newExps;
inputs.plotd.plotlevel='noplot';
    
inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = strcat('sim-',int2str(i));

% SIMULATION
inputs.ivpsol.ivpsolver='cvodes';
inputs.ivpsol.senssolver='cvodes';
inputs.ivpsol.rtol=1.0D-7;
inputs.ivpsol.atol=1.0D-7;
    
sim = AMIGO_SModel(inputs);

fid = fopen(resultFileName,'a');
fprintf(fid, '%f %f %f\n',IPTGe,sim.sim.states{1,1}(end,10),sim.sim.states{1,1}(end,10)/492); % IPTG and Citrine (AU)
fclose(fid); 
