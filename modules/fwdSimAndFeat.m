function [dur_mean,dur_sem,invCdf_mean,invCdf_sem,env_] = fwdSimAndFeat(mu,x,zeta,dt_infer,t_max,nAvg,...
                                                                   xPerc_dur,minBurstDuration,xPerc_cdf)
%%% 09-08-20    first revision
%%% Benoit Duchet, University of Oxford

%%% foward simulation of envelope models to compute average burst duration
%%% profiles and inverse cumulative distribution function.

%%% INPUTS
% mu:                   vector of drift function values
% x:                    vector of corresponding x values   
% zeta:                 noise standard deviation         
% dt_infer:             should be the value of dt_infer used for inference
% t_max:                duration of the data generated when forward simulating the model
% nAvg:                 number of repeats
% xPerc_dur:            vector of thresholds used for average burst
%                       duration profiles
% minBurstDuration:     should be set to 0. Bursts below this duration (in s) are not 
%                       considered when calculating average burst duration profiles. 
%                       The passage is expected to break down for values larger than 0.
% xPerc_cdf:            vector of thresholds used for inverse cdf

%%% OUTPUTS
% dur_mean:             average burst duration profile of foward simulated data
%                       (mean of nAvg repeats)
% dur_sem:              associated standard error of the mean
% invCdf_mean:          inverse cdf of foward simulated data
%                       (mean of nAvg repeats)
% invCdf_sem:           associated standard error of the mean
% env_:                 envelope forward simulated during the last repeat

x_min = x(1);
x_max = x(end);
dx = x(2) - x(1);
n = length(x);
theta_c = 17.6;
n_max = round(t_max/dt_infer);
nAvg_dur = 1;

for iAvg = 1:nAvg
    %%% envelope model forward simulation
    env(iAvg,:) = fwdSim_learntDyn_mex(x_min, x_max, dx, n, mu, zeta,...
                                       theta_c, n_max, dt_infer, normrnd(0,1,1,n_max-1));
            
    %%% average burst duration
    avgBurstDur(iAvg,:) = burstDurWrapper(env(iAvg,:),xPerc_dur,nAvg_dur,dt_infer,minBurstDuration,[]);
    
    %%% inverse cdf
    invCdf(iAvg,:) = prctile(env(iAvg,:),xPerc_cdf);
    
end  

%%% averaging
[invCdf_mean,invCdf_sem] = getMeanSem(invCdf,1);
[dur_mean,dur_sem] = getMeanSem(avgBurstDur,1);      
     
%%% returning last envelope segment
env_ = env(end,:);

end