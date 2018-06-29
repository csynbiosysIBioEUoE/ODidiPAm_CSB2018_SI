% This file requires the experimental data published in \ref.7 of the paper
% Create a structure with hysteresis data on yRG500 to be used for
% parameter estimation

Data = struct('expName',[],'t_con',[],'input',[],'n_samples',[],'t_samples',[],'exp_data',[],'iqr',[],'standard_dev',[]); 
Data.expName = {'up_0','up_2.5','up_5','up_7.5','up_10','up_15','up_20','up_25','up_35','up_50','up_100','up_1000',...
                'down_0','down_2.5','down_5','down_7.5','down_10','down_15','down_20','down_25','down_35','down_50','down_100','down_1000'};

for i=1:length(Data.expName)
    Name = strsplit(Data.expName{1,i},'_');

    if contains(Data.expName{1,i},'up')
        Data.input{1,i} = str2double(Name{1,2}); 
        Data.n_samples{1,i} = 5; 
        Data.t_samples{1,i} = [36 48 60 72 84].*3600; % sampling times in seconds
        Data.t_con{1,i} = [0 36 84].*3600; % tcon in seconds
        Data.exp_data{1,i} = normalised(find(normalised(:,1)==str2double(Name{1,2})),3);
        Data.iqr{1,i} = quan_75(find(quan_75(:,1)==str2double(Name{1,2})),3)-quan_25(find(quan_25(:,1)==str2double(Name{1,2})),3);
        Data.standard_dev{1,i} = Data.iqr{1,i}/(2*0.6745); 
    else
        Data.input{1,i} = str2double(Name{1,2}); 
        Data.n_samples{1,i} = 5; 
        Data.t_samples{1,i} = [36 48 60 72 84].*3600; % sampling times in seconds
        Data.t_con{1,i} = [0 36 84].*3600; % tcon in seconds
        Data.exp_data{1,i} = normalised_down(find(normalised_down(:,1)==str2double(Name{1,2})),3);
        Data.iqr{1,i} = quan_75d(find(quan_75d(:,1)==str2double(Name{1,2})),3)-quan_25d(find(quan_25d(:,1)==str2double(Name{1,2})),3);
        Data.standard_dev{1,i} = Data.iqr{1,i}/(2*0.6745); 
    end
end
%%
save('Ref7_Data.mat','Data')