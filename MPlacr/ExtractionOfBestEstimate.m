% Extract the iteration yielding the minimum SSE over the test set. 
Folder = '...'; % Specify the folder in which the mat data are located
cd   '...';     % Move to that folder

filePattern = fullfile(Folder, strcat('StellingFittingCrossValidationEqual','-','*.mat'));
Files = dir(filePattern); 
%%
SSE_Mat = zeros(length(Files),8);
for i=1:length(Files)
    FileName = Files(i).name;
    load(FileName);
    SSE_Mat(i,:) = SSE;
end
%%
SSE_vect = sum(SSE_Mat,2);
figure; 
plot([1:1:100]',SSE_vect);hold on; 
plot(find(SSE_vect == min(SSE_vect)),SSE_vect(find(SSE_vect == min(SSE_vect))),'*r')
xlabel('iteration index')
ylabel('SSE on the test set')
%% Load the file corresponding to the minimum SSE
load(Files(find(SSE_vect == min(SSE_vect))).name);
exps_indexTest = [16,23,2,18,1,12,5,9]
for i=1:sim_exps.n_exp
    figure;
    errorbar(sim_exps.t_s{1,i}/3600,sim_exps.exp_data{1,i},sim_exps.error_data{1,i},'ok'); hold on; 
    plot(sim_results.sim.tsim{1,i}/3600,sim_results.sim.obs{1,i},'b')
    title(int2str(exps_indexTest(i)))
    legend('experimental data', 'best fit')
    xlabel('Time (hours)')
    ylabel('Citrine (AU)')
end
