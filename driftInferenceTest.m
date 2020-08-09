%%% 09-08-20        first revision
%%% Benoit Duchet, University of Oxford

close all
clearvars
addpath(genpath(['.' filesep 'modules']))

%%% Parameters
toTest          = 'ou'; %choose from 'ou','deg2','deg3','deg4','deg5'
dt_infer        = 1E-3; %should be adapted to the roughness of the envelope
%                        use ~0.05 for the envelope of filtered LFP at beta frequency,
%                        0.001 for the synthetic test data provided (if using synthDat.(toTest).env directly)
pctToUse        = 50; %how much of the data given in env should be used for inference (give 50 for 50%)
plotDataMu      = true; %set to false if no groud truth is available for the data drift function (e.g. with patient data)
showActualZeta  = true; %same as above for the value of the noise standard deviation
runZetaOptim    = true; %run an optimisation to estimate the noise standard deviation,
%                        takes about two minutes. Running the script takes about three seconds without the optimisation
zeta_in         = []; %if above is set to false, zeta_in is used as the noise standard deviation

%%% loading data - should be adapted to use this script with your own datasets
temp = load('testMat');
synthDat = temp.synthDat;
SR = 1/synthDat.(toTest).dt; %%% this should be the sampling rate of your time series
timeSeries = synthDat.(toTest).z; %%% this should be your time series
env = synthDat.(toTest).env; %%% this should be the envelope of your time series
zeta = synthDat.(toTest).zeta; %%% remove if ground truth not available

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lower level parameters
w_cdf           = 0.1; %low weight given to the cdf in the optimisation cost. Only here to prevent
%                       the (rare) occurence of aberrant solutions.
n_x             = 300; %number of thresholds used by the passage method
minXfact        = 50; %determines the position of the lowest threshold used by the passage method as max(env) / minXfact
maxXfact        = 0.9; %determines the position of the highest threshold used by the passage method as max(env) * maxXfact
minBurstDuration = 0; %should be set to 0. Bursts below this duration (in s) are not considered by the passage
%                      method and when calculating average burst duration profiles. The passage is expected
%                      to break down for values larger than 0.
smoothMethod    = 'LOWESS'; %smoothing method used by the passage method. See the Matlab help
%                            of smoothdata for more details. Available methods are 'movmean',
%                            'movmedian', 'gaussian', 'lowess', 'loess', 'rlowess', 'rloess', 'sgolay'.
smoothFact1     = 5; %determines the smoothing span of T as floor( n_x / smoothFact1 )
smoothFact2     = 8; %determines the smoothing span of the derivative of T as floor( n_x / smoothFact2 )
interIdx        = [-1,-1];%determines how smoothing is handled at the left edge of T and its derivative.
%                          Negative values: keep all of the smoothed profiles, 0: manual selection (not recommanded if
%                          running the noise optimisation), 1E3: skips the lower 5% of the thresholds, then selects the
%                          first intersection, >0: idx of intersection to pick
xPerc_cdf       = 1:99; %vector of thresholds used in cdf plots and in the noise optimisation
xPerc_dur       = 20:5:95; %vector of thresholds used in average burst duration plots
maxFunCalls     = 40; %function call limit for the optimisation carried out to estimate the noise
zetaFact0_vect  = [0.05 0.2 0.5 0.8];%vector of initial noise values for local optimisations
n_muSem         = 3; %the data is divided in n_muSem to obtain sem error bars on the drift estimation
t_max           = 8E2; %duration of the data generated when forward simulating the inferred model
%                       to compare average burst duration profile and cdf to that of the original data
nAvg            = 5; %number of repeats for averaging when forward simulating the inferred model
showZetaPlot    = false; %to show the results of all local noise optimisations
mu_cutoff       = 98; %plot mu up to mu_cutoff percentile of x (note: 100% doesn't work)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% preparing data features and parameters for the noise optimisation
% only using the first pctToUse % of the data
N0 = length(env);
N1 = round(N0 * pctToUse / 100);
t1 = N1/SR;
env = env(1:N1);
timeSeries = timeSeries(1:N1);

% params object
params.x_min = max(env)/minXfact;
params.x_max = maxXfact*max(env);
params.n_x = n_x;
params.dt_infer = dt_infer;
params.minBurstDuration = minBurstDuration;
params.smoothMethod = smoothMethod;
params.smoothFact1 = smoothFact1;
params.smoothFact2 = smoothFact2;
params.interIdx = interIdx;
params.t_max = t_max;
params.nAvg = nAvg;
params.xPerc_dur = xPerc_dur;
params.xPerc_cdf = xPerc_cdf;
params.w_cdf = w_cdf;

%dat object
dat.env = env;
dat.SR = SR;
dat.zetaEff = std(timeSeries)/sqrt(1/SR)/3;
dat.dur_mean_dat = burstDurWrapper(dat.env,xPerc_dur,1,1/SR,minBurstDuration,[]); %average burst duration profile
dat.invCdf_mean_dat = prctile(dat.env,xPerc_cdf); %inverse cumulative distribution function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% optimisation to estimate the noise standard deviation (zeta)
if runZetaOptim
    singleArgCostFunc = @(x) costFunctionPassage(x, dat, params);
    options = optimset('MaxFunEvals',maxFunCalls);
   
    n_z = length(zetaFact0_vect);
    
    %n_z local optimisations
    for i_z = 1:n_z
        zetaFact0 = zetaFact0_vect(i_z);
        [z(i_z),fval(i_z),exitflag,output] = fminsearch(singleArgCostFunc,zetaFact0,options);
    end
    [~,i_min] = min(fval);
    bestZetaFact = z(i_min);
    res.zeta = bestZetaFact*dat.zetaEff;
    
    if showZetaPlot
        figure
        scatter(z*dat.zetaEff,fval)
        xlabel('\zeta')
        ylabel('cost')
    end
else
    res.zeta = zeta_in;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% passage method to get mu using the identified noise standard deviation
[res.mu,res.x] = passageMethod(env,SR,params.x_min,params.x_max,n_x,...
    res.zeta,dt_infer,minBurstDuration,...
    smoothMethod,smoothFact1,smoothFact2,interIdx);
%%% passage method on n_muSem segments to get error bars on mu
n1 = floor(length(env)/nAvg);
for i_k = 1:n_muSem
    env_k = env((1+n1*(i_k-1)):(n1*i_k));
    mu_k(i_k,:) = passageMethod(env_k,SR,params.x_min,params.x_max,n_x,...
        res.zeta,dt_infer,minBurstDuration,...
        smoothMethod,smoothFact1,smoothFact2,interIdx);
end
res.mu_sem = std(mu_k,1)/sqrt(n_muSem);

%%% obtaining model average burst duration profile and cdf
[res.dur_mean,res.dur_sem,res.invCdf_mean,res.invCdf_sem,res.env] = fwdSimAndFeat(...
    res.mu,res.x,res.zeta,dt_infer,t_max,nAvg,...
    xPerc_dur,minBurstDuration,xPerc_cdf);

%getting threshold indices up to mu_cutoff percentile of x
x_max_mod = res.invCdf_mean(find(xPerc_cdf>=mu_cutoff,1));
res.idx_mod = res.x < x_max_mod;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plotting the results
figure
subplot(2,4,[1 5])
hold on
if plotDataMu
    plot(synthDat.(toTest).x, synthDat.(toTest).mu,'color', 'k','LineWidth',1,'displayName','data')
end
plot(res.x(res.idx_mod), res.mu(res.idx_mod),'color', 'b','LineWidth',1,'displayName','inferred')
plotErr(res.x(res.idx_mod),res.mu(res.idx_mod)',res.mu_sem(res.idx_mod),'b',[],true)
legend('location','best')
xlabel('x')
ylabel('\mu (drift)')

title(['used ' num2str(t1) ' s'])

subplot(2,4,[2 6])
hold on
plot(xPerc_dur,dat.dur_mean_dat,'k','LineWidth',1,'displayName','data')
plot(xPerc_dur,res.dur_mean,'color','b','LineWidth',1,'displayName','inferred')
plotErr(xPerc_dur,res.dur_mean,res.dur_sem,'b',[],true)
xlabel('Threshold (%)')
ylabel('Average burst duration (s)')
legend('location','best')

if showActualZeta
    title(['actual \zeta = ' num2str(zeta,2)])
end

subplot(2,4,[3 7])
hold on
plot(xPerc_cdf,dat.invCdf_mean_dat,'k','LineWidth',1,'displayName','data')
plot(xPerc_cdf,res.invCdf_mean,'color','b','LineWidth',1,'displayName','inferred')
plotErr(xPerc_cdf,res.invCdf_mean,res.invCdf_sem,'b',[],true)
xlabel('Threshold (%)')
ylabel('x')
legend('location','best')

if runZetaOptim
    title(['inferred \zeta = ' num2str(res.zeta,2)])
else
    title(['\zeta = ' num2str(res.zeta,2)])
end
subplot(2,4,4)
n_dt_infer = 1000;
t_mod = (1:n_dt_infer)*dt_infer - dt_infer;
plot(t_mod,res.env(1:n_dt_infer),'b')
legend('inferred')
xlabel('Time (s)')
ylabel('x')

subplot(2,4,8)
n_dat = round(t_mod(end) * SR);
t_dat = (1:n_dat)/SR - 1/SR;
plot(t_dat,env(1:n_dat),'k')
legend('data')
xlabel('Time (s)')
ylabel('x')
