% This script simulates the response of MIPr to a prototype input taken
% from each class and extracts csv data used to create Fig. 3

cd ('../../');
AMIGO_Startup();

cd ('Examples/PLacCDC2018/MIPr/Scripts');

clear model;
clear exps;
clear best_global_theta;
clear pe_results;
clear pe_inputs;


% Specify folder name and short_name
results_folder = strcat('InduciblePromoter',datestr(now,'yyyy-mm-dd-HHMMSS'));
short_name     = strcat('IP');

% Load experimental data
load('PseudoData.mat');

% Read the model into the model variable
M3D_load_model; 

exps_indexall = [1,5,7,10];
exps.n_exp = length(exps_indexall);
inputs.model = model;


%Compute the steady state considering the initial theta guess and 0 IPTG
y0 = M3D_steady_state(inputs.model.par,0);

for iexp=1:length(exps_indexall)
    exp_indexData = exps_indexall(iexp);
    exps.exp_type{iexp} = 'fixed'; 
    exps.n_obs{iexp} = 1; 
    exps.obs_names{iexp} = char('Fluorescence');
    exps.obs{iexp} = char('Fluorescence = Cit_fluo');
    exps.u_interp{iexp} = Data.u_interp{1,exp_indexData};
    exps.t_f{iexp} = Data.t_con{1,exp_indexData}(end); 
    exps.n_s{iexp} = Data.n_s{1,exp_indexData};
    exps.t_s{iexp} = Data.t_s{1,exp_indexData}; 
    exps.t_con{iexp} = Data.t_con{1,exp_indexData};

    if exp_indexData>3 && exp_indexData<7 % Pulses
        exps.n_pulses{iexp} = Data.n_pulses{1,exp_indexData};
        exps.u_min{iexp} = Data.u_min{1,exp_indexData};
        exps.u_max{iexp} = Data.u_max{1,exp_indexData};
    else
        exps.n_steps{iexp} = length(exps.t_con{iexp})-1;
        exps.u{iexp} = Data.u{1,exp_indexData};
    end

    exps.data_type = 'real'; % check if this is the case
    exps.noise_type = 'homo_var'; % check if this is the case
    exps.exp_data{iexp} = Data.exp_data{1,exp_indexData};
    exps.error_data{iexp} = Data.error_data{1,exp_indexData};
    exps.exp_y0{iexp} = y0;
end

%% Compile and run the model 

inputs.exps  = exps;

inputs.pathd.results_folder = results_folder;
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       ='-simIPPseudo';
inputs.plotd.plotlevel='full';
% SIMULATION
inputs.ivpsol.ivpsolver='cvodes';
inputs.ivpsol.senssolver='cvodes';
inputs.ivpsol.rtol=1.0D-8;
inputs.ivpsol.atol=1.0D-8;
inputs.plotd.plotlevel='noplot';

AMIGO_Prep(inputs);
sim_results_Pseudo = AMIGO_SObs(inputs);


%% Now plotting the results

for e= 1:exps.n_exp
    figure; 
    errorbar(exps.t_s{1,e}(1:20:end)/60,exps.exp_data{1,e}(1:20:end),exps.error_data{1,e}(1:20:end),'ob'); hold on;
    plot(sim_results_Pseudo.sim.tsim{1,e}/60,sim_results_Pseudo.sim.obs{1,e},'b'); hold on;
    xlabel('Time (hours)');
    ylabel('Citrine (molecules)');
    title(int2str(e))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract data of input - output (simulation and pseudo-data)
% for each input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create csv step Data And Simulation
e = 1;
time_exp_step = exps.t_s{1,e}(1:20:end)/60; % hours
Citrine_molec_step = exps.exp_data{1,e}(1:20:end);
Citrine_molec_std_step = exps.error_data{1,e}(1:20:end);

T_exp = table(time_exp_step',Citrine_molec_step,Citrine_molec_std_step);
writetable(T_exp,'PseudoExperimentalData_Step.csv')
%
% Simulation
e = 1;
time_step = sim_results_Pseudo.sim.tsim{1,e}/60;
Citrine_molec_step = sim_results_Pseudo.sim.obs{1,e};
T_step = table(time_step',Citrine_molec_step);
writetable(T_step,'SimulationPseudo_Step.csv')
% Create csv pulse Data And Simulation
e = 2;
time_exp_pulse = exps.t_s{1,e}(1:20:end)/60; % hours
Citrine_molec_pulse = exps.exp_data{1,e}(1:20:end);
Citrine_molec_std_pulse = exps.error_data{1,e}(1:20:end);

T_exp = table(time_exp_pulse',Citrine_molec_pulse,Citrine_molec_std_pulse);
writetable(T_exp,'PseudoExperimentalData_Pulse.csv')
%
% Simulation
e = 2;
time_pulse = sim_results_Pseudo.sim.tsim{1,e}/60;
Citrine_molec_pulse = sim_results_Pseudo.sim.obs{1,e};
T_pulse = table(time_pulse',Citrine_molec_pulse);
writetable(T_pulse,'SimulationPseudo_Pulse.csv')
% Create csv pulse Data And Simulation
e = 3;
time_exp_ramp = exps.t_s{1,e}(1:20:end)/60; % hours
Citrine_molec_ramp = exps.exp_data{1,e}(1:20:end);
Citrine_molec_std_ramp = exps.error_data{1,e}(1:20:end);

T_exp = table(time_exp_ramp',Citrine_molec_ramp,Citrine_molec_std_ramp);
writetable(T_exp,'PseudoExperimentalData_Ramp.csv')
%
% Simulation
e = 3;
time_ramp = sim_results_Pseudo.sim.tsim{1,e}/60;
Citrine_molec_ramp = sim_results_Pseudo.sim.obs{1,e};
T_ramp = table(time_ramp',Citrine_molec_ramp);
writetable(T_ramp,'SimulationPseudo_Ramp.csv')
% Create csv pulse Data And Simulation
e = 4;
time_exp_random = exps.t_s{1,e}(1:20:end)/60; % hours
Citrine_molec_random = exps.exp_data{1,e}(1:20:end);
Citrine_molec_std_random = exps.error_data{1,e}(1:20:end);

T_exp = table(time_exp_random',Citrine_molec_random,Citrine_molec_std_random);
writetable(T_exp,'PseudoExperimentalData_Random.csv')

%
% Simulation
e = 4;
time_random = sim_results_Pseudo.sim.tsim{1,e}/60;
Citrine_molec_random = sim_results_Pseudo.sim.obs{1,e};
T_random = table(time_random',Citrine_molec_random);
writetable(T_random,'SimulationPseudo_Random.csv')

% 
figure; 
% input step 
IPTG_step = [repmat([zeros(1,250) 5*ones(1,250)],1,6)]; 
plot(1:1:3000,IPTG_step)

% input pulse
IPTG_pulse = [repmat([10*ones(1,50) 100*ones(1,10)],1,49) 10*ones(1,60)]; 
figure; 
plot(1:1:3000,IPTG_pulse)

% input ramp 
IPTG_ramp = [];
for i=1:round(length(exps.u{1,3})/2)
    IPTG_ramp = [IPTG_ramp exps.u{1,3}(i).*ones(1,60)];

end
IPTG_ramp_compl = [IPTG_ramp fliplr(IPTG_ramp(61:length(IPTG_ramp)-60))]

figure; 
plot(1:1:3000,IPTG_ramp_compl)

%%
IPTG_random = [];
for i=1:length(exps.u{1,4})-1
    IPTG_random = [IPTG_random exps.u{1,4}(i).*ones(1,150)];

end

figure; 
plot(1:1:3000,IPTG_random)

% Create csv with IPTG and time for each input

% Step 
IPTG_step = [repmat([zeros(1,250) 5*ones(1,250)],1,6)]; 
time_step = (1:1:3000)/60;

T_step = table(time_step',IPTG_step');
writetable(T_step,'Input_Step.csv')
%%
% Pulse 
IPTG_pulse = [repmat([10*ones(1,50) 100*ones(1,10)],1,49) 10*ones(1,60)]; 
time_pulse = (1:1:3000)/60;

T_pulse = table(time_pulse',IPTG_pulse');
writetable(T_pulse,'Input_Pulse.csv')
%%
% Ramp 
IPTG_ramp = IPTG_ramp_compl;
time_ramp = (1:1:3000)/60;

T_ramp = table(time_ramp',IPTG_ramp');
writetable(T_ramp,'Input_Ramp.csv')

%% Random 
time_random = (1:1:3000)/60;

T_random = table(time_random',IPTG_random');
writetable(T_random,'Input_Random.csv')