function plotErr(x,ymean,yerr,baseColor,alpha,noLeg)

%%% 26-11-19        checking the input vectors are rows
%%%                 option to skip legend
%%%                 alpha as an argument, and default value

%%% Benoit Duchet, University of Oxford

assert(isrow(x) && isrow(ymean) && isrow(yerr),'x, ymean and yerr have to be row vectors')
if isempty(alpha)
    alpha = 0.4;
end

idxToRm = isnan(ymean) | isnan(yerr);
x(idxToRm) = [];
ymean(idxToRm) = [];
yerr(idxToRm) = [];

if noLeg
    patch([x fliplr(x)],[ymean+yerr fliplr(ymean-yerr)],baseColor,'EdgeColor','none','FaceAlpha',alpha,'HandleVisibility','off');
else
   patch([x fliplr(x)],[ymean+yerr fliplr(ymean-yerr)],baseColor,'EdgeColor','none','FaceAlpha',alpha);
end 

end