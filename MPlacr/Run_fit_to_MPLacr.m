function [] = Run_fit_to_MPLacr( resultBase, numExperiments )
% This function runs the script for the identification of the model structure
% proposed in \ref 7 of the paper to their experimental data, using cross validation. It takes two inputs:
% - a string with the identifier of the output files
% - the number of iterations of parameter estimation. 
% In the paper we used numExperiments = 100.
cd ('../../');
AMIGO_Startup();

cd ('Examples/PLacCDC2018/MPLacr/Scripts');

% Specify the allowed boundaries for the parameters, the values were taken
% from \ref 7 of the paper
theta_min = [0.001,1e-10,7.75e-5 ,1e-10, 7.75e-5, 0.001, 5e-4, 1, 1,0.001, 1e-10,10];
theta_max = [1,1,1, 1, 1, 1, 1.7e-3, 1e4, 1000, 1, 1,100];

% Create a matrix of initial guesses for the parameters, having as many
% rows as the number of PE iterations (numExperiments)
% Each vector is passed as input to the computing function
M_norm = lhsdesign(numExperiments,length(theta_min));
M = zeros(size(M_norm));
for c=1:size(M_norm,2)
    for r=1:size(M_norm,1)
        %M(r,c) = M_norm(r,c)*(theta_max(1,c)-theta_min(1,c))+theta_min(1,c);
        M(r,c) = 10^(M_norm(r,c)*(log10(theta_max(1,c))-log10(theta_min(1,c)))+log10(theta_min(1,c))); % log exploration
    end
end 

ParFull = [M(:,1:2) 7.75e-5*ones(size(M,1),1) M(:,3:8) 2800*ones(size(M,1),1) M(:,9:end)];
save('MatrixParameters.mat','ParFull');

parfor epcc_exps=1:numExperiments
        try
            global_theta_guess = ParFull(epcc_exps,:);
            epccOutputResultFileNameBase = [resultBase,'-',num2str(epcc_exps)];
            [out] = fit_to_MPLacr(epccOutputResultFileNameBase,epcc_exps,global_theta_guess);
        catch err
            %open file
            errorFile = [resultBase,'-',num2str(epcc_exps),'.errorLog'];
            fid = fopen(errorFile,'a+');
            fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'));
            % close file
            fclose(fid);
        end
end

