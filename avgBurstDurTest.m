%%% 10-07-20        first revision
%%% Benoit Duchet, University of Oxford

%%% Example showing how to obtain the average burst duration profile of a
%%% time series using the code provided. This script is using the envelopes
%%% provided in testMat.

close all
clearvars
addpath(genpath(['.' filesep 'modules']))

%%% parameters
toTest              = 'deg5'; %choose from 'ou','deg2','deg3','deg4','deg5'
nAvg                = 5; % the time series is divided into nAvg segments to plot sem error bars 
minBurstDuration    = 0.1; % burst are only considered if longer than this duration (in s)
xPerc               = 20:5:95; % vector of thresholds

%%% loading data
temp = load('testMat');
synthDat = temp.synthDat;
dt = synthDat.(toTest).dt; %sampling interval
env = synthDat.(toTest).env; %envelope to analyse

%%% obtaining the average duration profile and errorbars.
[avgBurstDur,semBurstDur] = burstDurWrapper(env,xPerc,nAvg,dt,minBurstDuration,[]);
                                   
%%% plotting
figure
errorbar(xPerc,avgBurstDur,semBurstDur,'k--','displayName',toTest)
xlabel('Threshold (%)')
ylabel('Average burst duration (s)')
legend('location','best')
