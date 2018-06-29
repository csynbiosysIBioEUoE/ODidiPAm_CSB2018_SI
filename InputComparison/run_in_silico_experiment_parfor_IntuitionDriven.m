function [] = run_in_silico_experiment_parfor_IntuitionDriven(resultBase, numExperiments )

cd ('../../');
AMIGO_Startup();

cd ('Examples/PLacCDC2018/InputComparison/Scripts');

% Upper and lower boundaries selected for the parameters (theta refers to the parameter vector)
theta_min = [3.88e-5,3.88e-2,0.5,2,7.7e-3,0.2433,5.98e-5,0.012];
theta_max = [0.4950,0.4950,4.9,10,0.23,6.8067,0.2449,0.0217];

% Create a matrix of initial guesses for the vector of parameter estimates.
% The matrix has as many rows as the number of PE iterations
% (numExperiments).
% The vectors are obtained as latin hypercube samples within the boundaries selected above.
% Each vector is passed as input to the computing function

% M_norm = lhsdesign(numExperiments,length(theta_min));
% M = zeros(size(M_norm));
% for c=1:size(M_norm,2)
%     for r=1:size(M_norm,1)
%         M(r,c) = 10^(M_norm(r,c)*(log10(theta_max(1,c))-log10(theta_min(1,c)))+log10(theta_min(1,c))); % log exploration
%     end
% end 
% 
% ParFull = M; % in this case I am fitting all the values
% save('MatrixParameters_InputComparison.mat','ParFull');

% For the sake of a fair comparison among input classes, all simulations
% used the same matrix of initial theta guesses
load('MatrixParameters_InputComparison.mat');       

parfor epcc_exps=1:numExperiments
        try
            global_theta_guess = ParFull(epcc_exps,:);
            epccOutputResultFileNameBase = [resultBase,'-',num2str(epcc_exps)];
            % Uncomment the desired function for parameter estimation
            [out] = fit_to_InduciblePromoter_Step(epccOutputResultFileNameBase,epcc_exps,global_theta_guess); % How uses the homovar data (5%)
            %[out] = fit_to_InduciblePromoter_Pulse(epccOutputResultFileNameBase,epcc_exps,global_theta_guess); % How uses the homovar data (5%)
            %[out] = fit_to_InduciblePromoter_Random(epccOutputResultFileNameBase,epcc_exps,global_theta_guess); % How uses the homovar data (5%)

        catch err
            %open file
            errorFile = [resultBase,'-',num2str(epcc_exps),'.errorLog'];
            fid = fopen(errorFile,'a+');
            fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'));
            % close file
            fclose(fid);
        end
end
