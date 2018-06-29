% Extract the iteration yielding the minimum SSE over the test set. 
Folder = '../MIPr/Data/SimulationsMatFiles';
cd   '../MIPr/Data/SimulationsMatFiles';

filePattern = fullfile(Folder, strcat('Sim-InduciblePromoter_PseudoHomo','-','*.mat'));
Files = dir(filePattern); 
%% Create a matrix with the SSE values over the test set for each PE iteration
SSE_Mat = zeros(length(Files),4);
for i=1:length(Files)
    FileName = Files(i).name;
    load(FileName);
    SSE_Mat(i,:) = SSE;
end
%% Plot the values of the SSE, summed over the test set, for each iteration
SSE_vect = sum(SSE_Mat,2);
figure; 
plot([1:1:100]',SSE_vect);hold on; 
plot(find(SSE_vect == min(SSE_vect)),SSE_vect(find(SSE_vect == min(SSE_vect))),'*r')
xlabel('iteration index')
ylabel('SSE on the test set')
%% Load the file corresponding to the minimum SSE
load(Files(find(SSE_vect == min(SSE_vect))).name);
exps_indexTest = [4,9,10,12]
% Plot simulations over the test set for the best parameter estimate
for i=1:sim_exps.n_exp
    figure;
    errorbar(sim_exps.t_s{1,i}/60,sim_exps.exp_data{1,i},sim_exps.error_data{1,i},'ok'); hold on; 
    plot(sim_results.sim.tsim{1,i}/60,sim_results.sim.obs{1,i},'b','LineWidth',2)
    title(int2str(exps_indexTest(i)))
    legend('experimental data', 'best fit')
    xlabel('Time (hours)')
    ylabel('Citrine (mol)')
end
%%
disp('The best estimate is')
best_global_theta