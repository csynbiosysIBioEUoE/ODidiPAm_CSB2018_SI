%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In silico experiment OID script - runs PE, OED, mock experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out]=fit_to_InduciblePromoter_Optimised_valuesOnly_eSS(epccOutputResultFileNameBase,epccNumLoops,stepd,epcc_exps,global_theta_guess)

    resultFileName = [strcat(epccOutputResultFileNameBase),'.dat'];
    rng shuffle;
    rngToGetSeed = rng;

    % Write the header information of the .dat file in which the results of
    % PE (estimates, relative confidence intervals, residuals, relative
    % residuals and the time required for computation) will be stored.
    fid = fopen(resultFileName,'w');
    fprintf(fid,'HEADER DATE %s\n', datestr(datetime()));
    fprintf(fid,'HEADER RANDSEED %d\n',rngToGetSeed.Seed);
    fclose(fid);

    startTime = datenum(now);

    clear model;
    clear exps;
    clear best_global_theta_log;
    clear pe_results;
    clear ode_results;

    results_folder = strcat('InduciblePromoter',datestr(now,'yyyy-mm-dd-HHMMSS'));
    short_name     = strcat('IndProm',int2str(epcc_exps));

    % Read the model into the model variable
    InduciblePromoter_load_model;

    % We start with no experiments
    exps.n_exp=0;

    % Defining boundaries for the parameters (taken from scientific literature)
    global_theta_max = [0.4950,0.4950,4.9,10,0.23,6.8067,0.2449,0.0217];       % Maximum allowed values for the parameters
    global_theta_min = [3.88e-5,3.88e-2,0.5,2,7.7e-3,0.2433,5.98e-5,0.012];    % Minimum allowed values for the parameters
    global_theta_guess = global_theta_guess';


    % Specify the parameters to be calibrated.
    % The selection depends on the identifiability analysis preceding the
    % comparison: parameters that are not identifiable will fixed to the best
    % estimate available for them.
    % In our case, all parameters are identifiable.
    param_including_vector = [true,true,true,true,true,true,true,true];

    % Compile the model
    clear inputs;
    inputs.model = model;
    inputs.pathd.results_folder = results_folder;
    inputs.pathd.short_name     = short_name;
    inputs.pathd.runident       = 'initial_setup';
    AMIGO_Prep(inputs);

    % Fixed parts of the experiment
    totalDuration = 50*60;               % Duration in of the experiment (minutes)
    numLoops = epccNumLoops;             % Number of OID loops
    duration = totalDuration/numLoops;   % Duration of each loop (in our case the number is 1)
    stepDuration = stepd;                % Duration of each step (minutes). Note that the value passed to the function exceeds the response time, quantified in 80 mins for MPLac,r


for i=1:numLoops
    
    %Compute the steady state considering the initial theta guess and 0 IPTG
    y0 = InduciblePromoter_steady_state(global_theta_guess,0);
    
    
    % Need to determine the starting state of the next part of the
    % experiment we wish to design the input for. If there are multiple loops,
    % We get this state by simulating the model using the current best theta
    % for the duration for which we have designed input.
    if exps.n_exp == 0
        oid_y0 = [y0]; 
        best_global_theta = global_theta_guess;
    else
        % Simulate the experiment without noise to find end state
        clear inputs;
        inputs.model = model;
        inputs.model.par = best_global_theta;
        inputs.exps = exps;
        
        inputs.plotd.plotlevel='noplot';
        inputs.pathd.results_folder = results_folder;
        inputs.pathd.short_name     = short_name;
        inputs.pathd.runident       = strcat('sim-',int2str(i));
        
        sim = AMIGO_SData(inputs);
        
        oid_y0 = [sim.sim.states{1}(end,:)];
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimal experiment design
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear inputs;
    inputs.model = model;
    inputs.exps  = exps;
    format long g
    
    
    % Add new experiment that is to be designed
    inputs.exps.n_exp = inputs.exps.n_exp + 1;                      % Number of experiments
    iexp = inputs.exps.n_exp;                                       % Index of the experiment
    inputs.exps.exp_type{iexp}='od';                                % Specify the type of experiment: 'od' optimally designed
    inputs.exps.n_obs{iexp}=1;                                      % Number of observables in the experiment
    inputs.exps.obs_names{iexp}=char('Fluorescence');               % Name of the observables in the experiment
    inputs.exps.obs{iexp}=char('Fluorescence = Cit_fluo');          % Observation function
    
    
    % Fixed parts of the experiment
    inputs.exps.exp_y0{iexp}=oid_y0;                                % Initial conditions
    inputs.exps.t_f{iexp}=duration;                                 % Duration of the experiment (minutes)
    inputs.exps.n_s{iexp}=duration/5+1;                             % Number of sampling times - sample every 5 min
    
    % OED of the input
    inputs.exps.u_type{iexp}='od';
    inputs.exps.u_interp{iexp}='stepf';                             % Stimuli definition for experiment: 'stepf' steps of constant duration
    inputs.exps.n_steps{iexp}=round(duration/stepDuration);         % Number of steps in the input
    inputs.exps.t_con{iexp}=0:stepDuration:(duration);            % Switching times
    inputs.exps.u_min{iexp}=0*ones(1,inputs.exps.n_steps{iexp});    % Lower boundary for the input value
    inputs.exps.u_max{iexp}=1000*ones(1,inputs.exps.n_steps{iexp}); % Upper boundary for the input value
    
    inputs.PEsol.id_global_theta=model.par_names(param_including_vector,:);
    inputs.PEsol.global_theta_guess=transpose(global_theta_guess(param_including_vector));
    inputs.PEsol.global_theta_max=global_theta_max(param_including_vector);  % Maximum allowed values for the parameters
    inputs.PEsol.global_theta_min=global_theta_min(param_including_vector);  % Minimum allowed values for the parameters
    
    
    inputs.exps.noise_type='homo_var';           % Experimental noise type: Homoscedastic: 'homo'|'homo_var'(default)
    inputs.exps.std_dev{iexp}=[0.05];
    inputs.OEDsol.OEDcost_type='Dopt';
    
    
    % SIMULATION
    inputs.ivpsol.ivpsolver='cvodes';                     % [] IVP solver: 'cvodes'(default, C)|'ode15s' (default, MATLAB, sbml)|'ode113'|'ode45'
    inputs.ivpsol.senssolver='cvodes';                    % [] Sensitivities solver: 'cvodes'(default, C)| 'sensmat'(matlab)|'fdsens2'|'fdsens5'
    inputs.ivpsol.rtol=1.0D-8;                            % [] IVP solver integration tolerances
    inputs.ivpsol.atol=1.0D-8;
    
    % OPTIMIZATION
    %oidDuration=600;
    inputs.nlpsol.nlpsolver='eSS';
    inputs.nlpsol.eSS.maxeval = 5e4;
    inputs.nlpsol.eSS.maxtime = 30e3;
    inputs.nlpsol.eSS.local.solver = 'fminsearch'; 
    inputs.nlpsol.eSS.local.finish = 'fmincon';
    
    inputs.nlpsol.eSS.local.nl2sol.maxiter  =     500;     % max number of iteration
    inputs.nlpsol.eSS.local.nl2sol.maxfeval =     500;     % max number of function evaluation
    inputs.nlpsol.eSS.log_var=1:inputs.exps.n_steps{iexp};
    inputs.plotd.plotlevel='noplot';
    
    inputs.pathd.results_folder = results_folder;
    inputs.pathd.short_name     = short_name;
    inputs.pathd.runident       = strcat('oed-',int2str(i));
        
    oed_start = now;
    
    results = AMIGO_OED(inputs);
    oed_results{i} = results;
    oed_end = now;
    
    results.plotd.plotlevel = 'noplot';
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create a new experiment to simulate with the merged input
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear newExps;
    newExps.n_exp = 1;
    newExps.n_obs{1}=1;
    newExps.obs_names{1}=char('Fluorescence');
    newExps.obs{1}= char('Fluorescence = Cit_fluo');
    newExps.exp_y0{1}= [y0];
    
    newExps.t_f{1}=i*duration;
    newExps.n_s{1}=(i*duration)/5 + 1;
    newExps.t_s{1}=0:5:(i*duration);
    
    newExps.u_interp{1}='step';
    newExps.n_steps{1}=(i*duration)/stepDuration;
    newExps.t_con{1}=0:stepDuration:(i*duration);
    
    % Merge the input signal
    
    if exps.n_exp == 0
        newExps.u{1}=results.oed.u{results.oed.n_exp};
    else
        newExps.u{1}=[exps.u{1} results.oed.u{results.oed.n_exp}];
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mock an experiment
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear inputs;
    inputs.model = model;
    inputs.exps = newExps;
    
    inputs.exps.data_type='pseudo';
    inputs.exps.noise_type='homo_var';
    inputs.exps.std_dev{1}=[0.05];
    
    inputs.plotd.plotlevel='noplot';
    
    inputs.pathd.results_folder = results_folder;
    inputs.pathd.short_name     = short_name;
    inputs.pathd.runident       = strcat('sim-',int2str(epcc_exps),'-',int2str(i));
    
    sim = AMIGO_SData(inputs);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now we need to add this experiment output to newExps and copy
    % to exps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    newExps.exp_type{1}='fixed';
    if exps.n_exp == 0
        newExps.exp_data{1}   = sim.sim.exp_data{1};
        newExps.error_data{1} = sim.sim.error_data{1};
    else
        newExps.exp_data{1}   = [ exps.exp_data{1}   ; sim.sim.exp_data{1}((size(exps.exp_data{1},1)+1):end)];
        newExps.error_data{1} = [ exps.error_data{1} ; sim.sim.error_data{1}((size(exps.exp_data{1},1)+1):end)];
    end
    
    newExps.data_type='real';                                     % Type of data: 'pseudo'|'pseudo_pos'|'real'
    newExps.noise_type='homo_var';                                % Experimental noise type: Homoscedastic: 'homo'|'homo_var'(default)
    
    exps = newExps;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameter estimation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear inputs;
    inputs.model = model;
    inputs.exps  = exps;
    
    inputs.pathd.results_folder = results_folder;
    inputs.pathd.short_name     = short_name;
    inputs.pathd.runident       = strcat('pe-',int2str(epcc_exps),'-',int2str(i));
    
    % GLOBAL UNKNOWNS (SAME VALUE FOR ALL EXPERIMENTS)
    inputs.PEsol.id_global_theta=model.par_names(param_including_vector,:);
    inputs.PEsol.global_theta_guess=transpose(global_theta_guess(param_including_vector));
    inputs.PEsol.global_theta_max=global_theta_max(param_including_vector);  % Maximum allowed values for the paramters
    inputs.PEsol.global_theta_min=global_theta_min(param_including_vector);  % Minimum allowed values for the parameters
    
    % COST FUNCTION RELATED DATA
    inputs.PEsol.PEcost_type='lsq';        % 'lsq' (weighted least squares default) | 'llk' (log likelihood) | 'user_PEcost'
    inputs.PEsol.lsq_type='Q_expmax';
    
    % SIMULATION
    inputs.ivpsol.ivpsolver='cvodes';
    inputs.ivpsol.senssolver='cvodes';
    inputs.ivpsol.rtol=1.0D-7;
    inputs.ivpsol.atol=1.0D-7;
    
    
    % OPTIMIZATION
    inputs.nlpsol.nlpsolver='eSS';
    inputs.nlpsol.eSS.maxeval = 200000;
    inputs.nlpsol.eSS.maxtime = 5000;
    inputs.nlpsol.eSS.local.solver = 'lsqnonlin';  % nl2sol not yet installed on my mac
    inputs.nlpsol.eSS.local.finish = 'lsqnonlin';  % nl2sol not yet installed on my mac
    inputs.rid.conf_ntrials=500;
    
    inputs.plotd.plotlevel='noplot';
    
    pe_start = now;
    results = AMIGO_PE(inputs);
    pe_inputs{i} = inputs;
    pe_results{i} = results;
    pe_end= now;
    
    results.plotd.plotlevel = 'noplot';
    
    %Save the best theta
    best_global_theta(param_including_vector)=results.fit.thetabest;
    
    %Write some results to the output file
    fid = fopen(resultFileName,'a');
    used_par_names = model.par_names(param_including_vector,:);
    
    for j=1:size(used_par_names,1)
        fprintf(fid,'ITERATION %d PARAM_FIT %s %f\n', i, used_par_names(j,:), results.fit.thetabest(j));
        if isfield(results.fit,'rel_conf_interval')
            fprintf(fid,'ITERATION %d REL_CONF %s %f\n',  i, used_par_names(j,:), results.fit.rel_conf_interval(j));
        end
        if isfield(results.fit,'residuals')
            fprintf(fid,'ITERATION %d RESIDUAL %s %f\n', i, used_par_names(j,:), results.fit.residuals{1}(j));
        end
        if isfield(results.fit,'rel_residuals')
            fprintf(fid,'ITERATION %d REL_RESIDUAL %s %f\n', i, used_par_names(j,:), results.fit.rel_residuals{1}(j));
        end
    end
    %Time in seconds
    fprintf(fid,'ITERATION %d OED_TIME %.1f\n', i, (oed_end-oed_start)*24*60*60);
    fprintf(fid,'ITERATION %d PE_TIME %.1f\n',  i, (pe_end-pe_start)*24*60*60);
    fclose(fid);

    
end

% Now log stuff at the end
for i=1:10

    duration = i*5*60;  % Duration in minutes

    clear inputs;
    inputs.model = model;
    inputs.exps  = exps;

    % Reduce the input to a smaller set of values
    inputs.exps.t_f{1}          = duration;                % Experiment duration
    inputs.exps.n_s{1}          = duration/5 + 1;          % Number of sampling times
    inputs.exps.t_s{1}          = 0:5:duration;            % Times of samples

    inputs.exps.n_steps{1}      = sum(exps.t_con{1} < duration);
    inputs.exps.t_con{1}        = exps.t_con{1}(1:inputs.exps.n_steps{1}+1);
    inputs.exps.u{1}            = exps.u{1}(1:inputs.exps.n_steps{1});

    inputs.exps.exp_data{1}     = exps.exp_data{1}(1:inputs.exps.n_s{1});
    inputs.exps.error_data{1}   = exps.error_data{1}(1:inputs.exps.n_s{1});

    inputs.pathd.results_folder = results_folder;
    inputs.pathd.short_name     = short_name;
    inputs.pathd.runident       = strcat('pe-',int2str(i));

    % GLOBAL UNKNOWNS (SAME VALUE FOR ALL EXPERIMENTS)
    inputs.PEsol.id_global_theta=model.par_names(param_including_vector,:);
    inputs.PEsol.global_theta_guess=transpose(global_theta_guess(param_including_vector));
    inputs.PEsol.global_theta_max=global_theta_max(param_including_vector);  % Maximum allowed values for the paramters
    inputs.PEsol.global_theta_min=global_theta_min(param_including_vector);  % Minimum allowed values for the parameters

    % COST FUNCTION RELATED DATA
    inputs.PEsol.PEcost_type='lsq';        % 'lsq' (weighted least squares default) | 'llk' (log likelihood) | 'user_PEcost'
    inputs.PEsol.lsq_type='Q_I';

    % SIMULATION
    inputs.ivpsol.ivpsolver='cvodes';
    inputs.ivpsol.senssolver='cvodes';
    inputs.ivpsol.rtol=1.0D-8;
    inputs.ivpsol.atol=1.0D-8;


    % OPTIMIZATION
    inputs.nlpsol.nlpsolver='eSS';
    inputs.nlpsol.eSS.maxeval = 200000;
    inputs.nlpsol.eSS.maxtime = 5000;
    inputs.nlpsol.eSS.local.solver = 'lsqnonlin';  % nl2sol not yet installed on my mac
    inputs.nlpsol.eSS.local.finish = 'lsqnonlin';  % nl2sol not yet installed on my mac
    inputs.rid.conf_ntrials=500;

    inputs.plotd.plotlevel='noplot';


    pe_start = now;
    results = AMIGO_PE(inputs);
    pe_inputs2{i} = inputs;
    pe_results2{i} = results;
    pe_end= now;

    % Write some results to the output file
    fid = fopen(resultFileName,'a');
    used_par_names = model.par_names(param_including_vector,:);

    for j=1:size(used_par_names,1)
        fprintf(fid,'HOUR %d PARAM_FIT %s %f\n', i*5, used_par_names(j,:), results.fit.thetabest(j));
        if isfield(results.fit,'rel_conf_interval')
            fprintf(fid,'HOUR %d REL_CONF %s %f\n',  i*5, used_par_names(j,:), results.fit.rel_conf_interval(j));
        end
        if isfield(results.fit,'residuals')
           fprintf(fid,'HOUR %d RESIDUAL %s %f\n', i*5, used_par_names(j,:), results.fit.residuals{1}(j));
        end
        if isfield(results.fit,'rel_residuals')
            fprintf(fid,'HOUR %d REL_RESIDUAL %s %f\n', i*5, used_par_names(j,:), results.fit.rel_residuals{1}(j));
        end
    end
    % Time in seconds
    fprintf(fid,'HOUR %d PE_TIME %.1f\n',  i*5, (pe_end-pe_start)*24*60*60);
    fclose(fid);

    best_global_theta_log{i}=results.fit.thetabest;

end

true_param_values = model.par(param_including_vector);

save([strcat(epccOutputResultFileNameBase),'.mat'], 'pe_inputs','pe_results','pe_results2','oed_results','exps','inputs','true_param_values','best_global_theta','best_global_theta_log');

out= 1;
end
