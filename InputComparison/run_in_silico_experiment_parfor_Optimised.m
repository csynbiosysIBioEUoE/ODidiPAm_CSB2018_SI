function [] = run_in_silico_experiment_parfor_Optimised( resultBase, numLoops, numExperiments )

cd ('../../');
AMIGO_Startup();

cd ('Examples/In_Silico_Loop');

% Selected boundaries for the parameters
theta_min = [3.88e-5,3.88e-2,0.5,2,7.7e-3,0.2433,5.98e-5,0.012];
theta_max = [0.4950,0.4950,4.9,10,0.23,6.8067,0.2449,0.0217];

% Create a matrix of initial guesses for the parameters, having as many
% rows as the number of PE iterations (numExperiments) 
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

load('MatrixParameters_InputComparison.mat');   

parfor epcc_exps=1:numExperiments
        stepd = 200;
        epccNumLoops = numLoops;
        try
            global_theta_guess = ParFull(epcc_exps,:);
            epccOutputResultFileNameBase = [resultBase,'-','OptstepseSS-',num2str(numLoops),'_loops-',num2str(epcc_exps)];
            [out]=fit_to_InduciblePromoter_Optimised_valuesOnly_eSS(epccOutputResultFileNameBase,epccNumLoops,stepd,epcc_exps,global_theta_guess);

        catch err
            %open file
            errorFile = [resultBase,'-','OptstepseSS-',num2str(numLoops),'_loops-',num2str(epcc_exps),'.errorLog'];
            fid = fopen(errorFile,'a+');
            fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'));
            % close file
            fclose(fid);
        end
end

