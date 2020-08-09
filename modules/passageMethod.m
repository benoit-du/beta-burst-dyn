function [mu,xTprime] = passageMethod(env,SR,...
                                      x_min,x_max,n_x,...
                                      zeta,dt_infer,minBurstDuration,...
                                      smoothMethod,smoothFact1,smoothFact2,interIdx)
                                  
%%% 09-08-20    first revision
%%% Benoit Duchet, University of Oxford

%%% implementation of the passage method to infer the envelope model drift
%%% function of a given envelope.

%%% INPUTS
% env:                  data envelope             
% SR:                   data sampling rate
% x_min:                lowest threhold to use
% x_max:                highest threshold to use
% n_x:                  number of thresholds to use
% zeta:                 noise standard deviation
% dt_infer:             should be adapted to the roughness of the envelope.
%                       Use ~0.05 for the envelope of filtered LFP at beta frequency,
%                       0.001 for the synthetic test data provided (if using synthDat.(toTest).env directly)
% minBurstDuration:     should be set to 0. Bursts below this duration (in s) are not considered. 
%                       The passage method is expected to break down for values larger than 0.
% smoothMethod:         smoothing method used by the passage method. See the Matlab help
%                       of smoothdata for more details. Available methods are 'movmean',
%                       'movmedian', 'gaussian', 'lowess', 'loess', 'rlowess', 'rloess', 'sgolay'.
% smoothFact1:          determines the smoothing span of T as floor( n_x / smoothFact1 )
% smoothFact2:          determines the smoothing span of T as floor( n_x / smoothFact2 )
% interIdx:             1x2 array which determines how smoothing is handled at the left edge of T and its derivative.
%                       Negative values: keep all of the smoothed profiles, 0: manual selection (not recommanded if
%                       running the noise optimisation), 1E3: skips the lower 5% of the thresholds, then selects the
%                       first intersection, >0: idx of intersection to pick

%%% OUTPUTS
% mu:                   inferred drift function    
% xTprime:              corresponding x vector


%%% obtaining average burst duration profiles and T
x = linspace(x_min,x_max,n_x); %learning range
avgBurstDuration = NaN(n_x,1);

for p = 1:n_x
    tempDur = getBurstDurationAmplitude_mex(env,length(env),SR,x(p),minBurstDuration);
    avgBurstDuration(p) = mean(tempDur(tempDur~=0));
end

T = sqrt(2/(pi*dt_infer)) / zeta * avgBurstDuration;

%%% smoothing T
smoothSp1 = floor( n_x / smoothFact1 );
smoothSp2 = floor( n_x / smoothFact2 );
T_sm = smoothdata(T,smoothMethod,smoothSp1);

%%% handling smoothing at the left edge of T based on interIdx
[x0,y0,iout] = intersections(x,T,x,T_sm);
if isempty(interIdx)
    interIdx = [1,1];
elseif interIdx(1) == 1E3
    interIdx(1) = find(x0>max(x)*0.05,1);
elseif interIdx(1) == 0
    figure
    hold on
    plot(x,T,'displayName','raw')
    plot(x,T_sm,'displayName','smoothed')
    for i_split = 1:length(iout)
        scatter(x0(i_split),y0(i_split),'handleVisibility','off')
    end
    xlabel('x','interpreter','latex')
    ylabel('T','interpreter','latex')
    legend('location','best')
    title('select intersection')
    [x_in,~] = ginput(1);
    close
    if x_in < 0
        i_selected = 0;
    else
        i_selected = getClosestIdx(x_in,x0);
    end
    interIdx(1) = i_selected;
end

if ~isempty(iout) && interIdx(1)>0
    n_split = floor(iout(interIdx(1)));
    T_sm(1:n_split) = T(1:n_split);
end

%%% derivative of T and smoothing
Tprime_fromSm = diff(T_sm)/(x(2)-x(1));
xTprime = x(1:end-1);
T_sm_ = T_sm(1:end-1);
Tprime_fromSm_sm = smoothdata(Tprime_fromSm,smoothMethod,smoothSp2);

%%% handling smoothing at the left edge of the derivative of T based on interIdx
[x0prime,y0prime,ioutprime] = intersections(xTprime,Tprime_fromSm,xTprime,Tprime_fromSm_sm);
if interIdx(2) == 1E3
    interIdx(2) = find(x0prime>max(xTprime)*0.05,1);
elseif interIdx(2) == 0
    figure
    hold on
    plot(xTprime,Tprime_fromSm,'displayName','raw')
    plot(xTprime,Tprime_fromSm_sm,'displayName','smoothed')
    for i_split = 1:length(ioutprime)
        scatter(x0prime(i_split),y0prime(i_split),'handleVisibility','off')
    end
    xlabel('x','interpreter','latex')
    ylabel('$T^\prime$','interpreter','latex')
    legend('location','best')
    title('select intersection')
    [x_in,~] = ginput(1);
    close
    if x_in < 0
        i_selected = 0;
    else
        i_selected = getClosestIdx(x_in,x0prime);
    end
    interIdx(2) = i_selected;
end

if ~isempty(ioutprime) && interIdx(2)>0
    n_splitPrime = floor(ioutprime(interIdx(2)));
    Tprime_fromSm_sm(1:n_splitPrime) = Tprime_fromSm(1:n_splitPrime);
end

%%% drift function
mu  = - ( 1 + 0.5 * zeta^2 * Tprime_fromSm_sm) ./ T_sm_;

end
