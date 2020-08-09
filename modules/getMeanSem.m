function [yMean,ySem] = getMeanSem(y,dim)

%%% 09-08-20    first revision
%%% Benoit Duchet, University of Oxford

%%% returning the mean and sem of y along dimension dim. Robust to NaNs.

yMean = nanmean(y,dim);
nObs = sum(~isnan(y),dim);
nObs(nObs == 0) = 1;
ySem = nanstd(y,[],dim)./sqrt(nObs);

end