function [] = Run_MPLac_in_silico_experiment( resultBase, numExperiments )
% This function calls MPLac_in_silico.m to simulate the MPLac model to a
% step in IPTG (the value of which is specified by IPTGExt) lasting 24
% hours. It has been used to simulate the MPLac dose response curve as well
% as the response of the system to a step of IPTG with amplitude 5uM. 
global epccOutputResultFileNameBase;
global IPTGe

cd ('../../');
AMIGO_Startup();

cd ('Examples/PLacCDC2018/MPLac/Scripts');
% a) Simulation of the dose-response curve
IPTGExt = [0:0.1:49, 50:10:1000,1100:1000:1e4];    

% b) Simulation of the dynamic to a step to 5 uM of IPTG (Figure 2b)
%IPTGExt = 5;
for epcc_exps=1:numExperiments
    for ip = 1:length(IPTGExt)
        IPTGe = IPTGExt(ip);
        try
            epccOutputResultFileNameBase = [resultBase,'-',num2str(epcc_exps)];
            MPLac_in_silico;
        catch err
            %open file
            errorFile = [resultBase,'-',num2str(epcc_exps),'.errorLog'];
            fid = fopen(errorFile,'a+');
            fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'));
            % close file
            fclose(fid);
        end
    end
end

