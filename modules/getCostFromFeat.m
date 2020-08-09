function cost = getCostFromFeat(model,data)

%%% 09-08-20    first revision
%%% Benoit Duchet, University of Oxford

%%% calculating the cost (1 - R squared) associated with one (vector) feature.

%%% INPUTS
% model:            vector value of the feature evaluted on the model output     
% data:             vector value of the feature evaluted on the data output                  

model = model(:);
data = data(:);

cost = sum((data-model).^2)/sum((data-mean(data)).^2);

end