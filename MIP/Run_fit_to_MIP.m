function [] = Run_fit_to_MIP( resultBase, numExperiments )
% This function runs the script for the identification of the M3D structure
% to the Stelling experimental data (\ref 7 of the paper). It takes two inputs:
% - a string with the identifier of the output files
% - the number of iterations of parameter estimation. 
% In the paper we used numExperiments = 100.

cd ('../../');
AMIGO_Startup();

cd ('Examples/PLacCDC2018/MIP/Scripts');

% Specify the allowed boundaries for the parameters (parameters expressed in 1/s as. The last parameter is taken from ref 7 of the paper)
theta_min = [6.4667e-7,6.4667e-4,0.5,2,1.2833e-4,0.0041,9.9667e-7,2e-4,10]; %sec
theta_max = [0.0083,0.0083,4.9,10,0.0038,0.1134,0.0041,3.6167e-4,100];


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
% save('MatrixParameters.mat','ParFull');

load('MatrixParameters.mat');   

parfor epcc_exps=1:numExperiments
        try
            global_theta_guess = ParFull(epcc_exps,:);
            epccOutputResultFileNameBase = [resultBase,'-',num2str(epcc_exps)];
            [out] = fit_to_MIP(epccOutputResultFileNameBase,epcc_exps,global_theta_guess);
        catch err
            %open file
            errorFile = [resultBase,'-',num2str(epcc_exps),'.errorLog'];
            fid = fopen(errorFile,'a+');
            fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'));
            % close file
            fclose(fid);
        end
end

