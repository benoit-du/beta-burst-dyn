function [avgBurstDuration,semBurstDuration,segAvgBurstDuration,perc] = burstDurWrapper(env,xPerc,nAvg,dt,minBurstDuration,percIn)

%%% 09-08-20        first revision
%%% Benoit Duchet, University of Oxford

%%% Mex function wrapper to obtain the average burst duration profile of an envelope

%%% INPUTS
% env:                 envelope time series
% xPerc:               vector of thresholds (format >> 75% should be given as 75). Can be empty if percIn is non empty.
% nAvg:                the time series is divided into nAvg segments to obtain sem error bars. Set to 1 if error bars are not needed.
% dt:                  time series sampling interval (should be in the same units as minBurstDuration).
% minBurstDuration:    bursts are not considered if shorter than minBurstDuration (should be in the same units as dt).
% percIn:              vector of thresholds (as envelope values rather than %). Should be empty if using xPerc.

%%% OUTPUTS
% avgBurstDuration:    average burst duration profile
% semBurstDuration:    sem error bars
% burstDuration:       matrix of average burst duration - rows correspond to segments, columns to thresholds
% perc:                vector of thresholds used (as envelope values rather than %)

if isempty(percIn)
    perc = prctile(env,xPerc);
else
    perc = percIn;
end

n1 = floor(length(env)/nAvg);

for k = 1:nAvg
    env_k = env((1+n1*(k-1)):(n1*k));
    for p = 1:length(perc)
        % a vector of duration of individual bursts longer than minBurstDuration is returned in tempDur
        % the mex function also returns a vector of peak amplitude of individual bursts (second output, not used here)
        tempDur = getBurstDurationAmplitude_mex(env_k,length(env_k),1/dt,perc(p),minBurstDuration);
        segAvgBurstDuration(k,p) = mean(tempDur(tempDur~=0));
    end
    
    
end

if k>1
    avgBurstDuration = nanmean(segAvgBurstDuration);
    semBurstDuration = nanstd(segAvgBurstDuration)/sqrt(size(segAvgBurstDuration,1));
else
    avgBurstDuration = segAvgBurstDuration;
    semBurstDuration = NaN(size(segAvgBurstDuration));
end
end