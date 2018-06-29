%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulations to compute SSE for the MPLacr model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write the header information
startTime = datenum(now);

clear model;
clear exps;
clear best_global_theta;
clear pe_results;
clear pe_inputs;

% Specify folder name and short_name
results_folder = strcat('PLac_OpenLoop',datestr(now,'yyyy-mm-dd-HHMMSS'));
short_name     = strcat('PLacOL');

% Load experimental data. This matfile was derived from the data by Gnugge
% et al., using the script Ref7_Data_StructureCreation.m. Please refer to
% that publication to access the original experimental data. 
load('Ref7_Data.mat');

% Read the model into the model variabl
PLac_load_model; 


%Compute the steady state considering the initial theta guess and 0 IPTG
y0 = PLac_Compute_SteadyState(model.par,0);

exps_indexall = [22,3,10,19,17,15,4,14,6,8,21,20,13,24,7,11,16,23,2,18,1,12,5,9];

exps.n_exp = length(exps_indexall);

for iexp=1:length(exps_indexall)

    exp_indexData = exps_indexall(iexp);
    exps.exp_type{iexp} = 'fixed'; 
    exps.n_obs{iexp} = 1; 
    exps.obs_names{iexp} = char('Fluorescence');
    exps.obs{iexp} = char('Fluorescence = Cit_AU');
    y0_to_use = PLac_Compute_SteadyState(model.par,0);   
    exps.exp_y0{iexp}=y0_to_use;  
    exps.t_f{iexp} = Data.t_con{1,exp_indexData}(end); 
    exps.n_s{iexp} = Data.n_samples{1,exp_indexData};
    exps.t_s{iexp} = Data.t_samples{1,exp_indexData}; 
    exps.u_interp{iexp} = 'step';
    exps.t_con{iexp} = Data.t_con{1,exp_indexData};
    exps.n_steps{iexp} = length(exps.t_con{iexp})-1;

    if exp_indexData<13
       exps.u{iexp} =  [0 Data.input{1,exp_indexData}];

    else
       exps.u{iexp} = [1000 Data.input{1,exp_indexData}];
    end

    exps.exp_data{iexp} = Data.exp_data{1,exp_indexData};
    exps.error_data{iexp} = Data.standard_dev{1,exp_indexData};
end

inputs.model = model;
inputs.exps  = exps;

inputs.pathd.results_folder = results_folder;
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       ='-sim';
inputs.plotd.plotlevel='noplot';

AMIGO_Prep(inputs);

% SIMULATION
inputs.ivpsol.ivpsolver='cvodes';
inputs.ivpsol.senssolver='cvodes';
inputs.ivpsol.rtol=1.0D-8;
inputs.ivpsol.atol=1.0D-8;

%inputs.plotd.plotlevel='full';
sim_results = AMIGO_SObs(inputs);

% Now compute SSE
SSE = zeros(size(exps_indexall));
for iexp=1:length(exps_indexall)
    exp_indexData = exps_indexall(iexp);
    SSE(iexp) = sum((Data.exp_data{1,exp_indexData}-sim_results.sim.sim_data{1,iexp}).^2);
end


sim_inputs = inputs;
sim_exps = exps;

save('ResultsSSE_StellingCrossValidationEqual.mat','sim_results','sim_inputs','sim_exps','SSE');


