function cost = costFunctionPassage(zetaFact, dat, params)

%%% 09-08-20    first revision
%%% Benoit Duchet, University of Oxford

%%% cost function used as part of the optimisation to estimate the noise standard
%%% deviation.

%%% INPUTS
% zetaFact:            value of the zeta factor to estimate the cost of. zetaFact 
%                      is related to zeta by zeta = zetaEff * zetaFact        
% dat:                 data object (see required fields in code below)
% params:              parameter object (see required fields in code below)
%
%%% OUTPUTS
% cost:                cost value

%%% unpacking the parameter object
x_min = params.x_min;
x_max = params.x_max;
n_x = params.n_x;
dt_infer = params.dt_infer;
minBurstDuration = params.minBurstDuration;
smoothMethod = params.smoothMethod;
smoothFact1 = params.smoothFact1;
smoothFact2 = params.smoothFact2;
interIdx = params.interIdx;
t_max = params.t_max;
nAvg = params.nAvg;
xPerc_dur = params.xPerc_dur;
xPerc_cdf = params.xPerc_cdf;
w_cdf = params.w_cdf;

%%% unpacking the data object
env = dat.env;
SR = dat.SR;
zetaEff = dat.zetaEff;
dur_mean_dat = dat.dur_mean_dat;
invCdf_mean_dat = dat.invCdf_mean_dat;

% obtaining zeta from zetaFact
zeta = zetaEff * zetaFact;

%%% passage method
[mu,x] = passageMethod(env,SR,...
                       x_min,x_max,n_x,...
                       zeta,dt_infer,minBurstDuration,...
                       smoothMethod,smoothFact1,smoothFact2,interIdx);
                   
%%% model features
[dur_mean,~,invCdf_mean] = fwdSimAndFeat(mu,x,zeta,dt_infer,t_max,nAvg,...
                                         xPerc_dur,minBurstDuration,xPerc_cdf);

%%% calculating the associated cost                                 
cost = (1-w_cdf)*getCostFromFeat(dur_mean,dur_mean_dat) + w_cdf*getCostFromFeat(invCdf_mean,invCdf_mean_dat);

end