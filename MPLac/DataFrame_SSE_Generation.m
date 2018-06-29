% Script to create and export the structure for SSE to be used in python
% notebook 
clc, clear all, close all;
% import data (the files have been generated running ComputeESSStellingVsData.m)
% Results of SSE using the Stelling parameters (MPLac)
load('ResultsSSE_StellingBestOL.mat');
SSE_Stelling = SSE';
% Results of SSE using the Stelling model, our fitted parameters (MPLac,r)
load('ResultsSSE_StellingCrossValidationEqual.mat');
SSE_Us = SSE';
% Results of SSE using the Inducible Promoter model, fitted on the stelling
% data (MIP)
load('ResultsSSE_StellingInduciblePromoter.mat');
SSE_IP = SSE';
% Create csv file
exptIndex_set = [22,3,10,19,17,15,4,14,6,8,21,20,13,24,7,11,16,23,2,18,1,12,5,9]';
type = {};
type = [repelem({'Stelling'},length(exptIndex_set))'];
type = [type; repelem({'Us'},length(exptIndex_set))'];
type = [type; repelem({'IP'},length(exptIndex_set))'];
exptIndex = repmat(exptIndex_set,3,1);
SSE = [SSE_Stelling; SSE_Us; SSE_IP];

index = linspace(1,length(exptIndex),length(exptIndex));
rowsi = strread(num2str(index),'%s');

DataSSE = table(type, exptIndex,SSE,'RowNames',rowsi);
writetable(DataSSE,'DataSSE.csv')