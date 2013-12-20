function [allSses, fit_complexity_vec, fit_criteria] = metadata_fit_models_crossval(dataInputStruct)

responseVars = dataInputStruct.responseVars;
dataUseInds = dataInputStruct.dataUseInds;
predictorVarsCat = dataInputStruct.predictorVarsCat;
predictorVarsCont = dataInputStruct.predictorVarsCont;
% predictorVars = dataInputStruct.predictorVars;
predictorVarsNames = dataInputStruct.predictorVarsNames;
responseVarsNames = dataInputStruct.responseVarsNames;
catVarInds = dataInputStruct.catVarInds;
NeuronType = dataInputStruct.NeuronType;

MAXFITSTEPS = 250;
fit_complexity_vec = {'constant','linear', 'purequadratic', 'interactions','quadratic'};
fit_criteria = {'bic', 'aic'};
numResamples = 100;
removeFract = .1;

% responseVars = [rmp	log10(ir) log10(tau)	amp	hw thresh];
% dataUseInds = log10(Age) > .67 & ~strcmp('cell culture ', PrepType);
% predictorVarsCat = [NeuronType	Species	Strain	ElectrodeType PrepType JxnPotential];
% predictorVarsCont = [log10(Age) Temp];

% predictorVarsNames = {'NeuronType';	'Species';	'Strain';	'ElectrodeType'; 'PrepType'; 'JxnPotential'; 'Age'; 'Temp'};
% responseVarsNames = {'rmp'; 'ir'; 'tau'; 'amp'; 'hw'; 'thresh'};

% catVarInds = [1 1 1 1 1 1 0 0];

predictorVars = [];
for i = 1:size(predictorVarsCat, 2)
    predictorVars(:,end+1) = toCategorical(predictorVarsCat(dataUseInds,i));
end

for i = 1:size(predictorVarsCont, 2)
    predictorVars(:,end+1) = predictorVarsCont(dataUseInds, i);
end

numResponseVars = size(responseVars, 2);
[neuronCategorical, neuronNames] = grp2idx(NeuronType(dataUseInds));
numUniqueNeurons = max(neuronCategorical);
minNeuronMentions = 1;
allSses = zeros(6, length(fit_complexity_vec), length(fit_criteria), numResamples);

for i = [1 2 3 4 5 6]
    for j = 1:length(fit_complexity_vec)
        complexity = fit_complexity_vec{j};
        for k = 1:length(fit_criteria)
            criteria = fit_criteria{k};
            parfor n = 1:numResamples
                
                %     if i == 3
                %         criteria = 'bic';
                %         complexity = 'linear';
                % %         complexity = fit_complexity;
                %     else
                %         criteria = fit_criteria;
                %         complexity = fit_complexity;
                %     end
                
                
                ephysInd = i;
                currRespVar = responseVars(dataUseInds,ephysInd);
                currUsePredVarInds = 1:length(predictorVarsNames);
                
                
                
                % only use data from neuron types which have at least 3 reports
                predRespMat = [predictorVars(:,currUsePredVarInds) currRespVar];
                predRespVec = sum(predRespMat, 2);
                articleInds = [];
                numUsedNeurons = 0;
                for m = 1:numUniqueNeurons
                    %     numUniqueArts = sum(neuronCategorical == m & ~isnan(currRespVar));
                    numUniqueArts = sum(neuronCategorical == m );
                    if numUniqueArts >= minNeuronMentions
                        %                     if randn > 0
                        %             currArticleInds = find(neuronCategorical == m & ~isnan(currRespVar));
                        currArticleInds = find(neuronCategorical == m);
                        articleInds = [articleInds; currArticleInds];
                        numUsedNeurons = numUsedNeurons + 1;
                        %                     end
                    end
                end
                
                notNaNInds = find(~isnan(currRespVar(articleInds)));
                randRemoveInds = randsample(notNaNInds, floor(removeFract*length(notNaNInds)));
                currRespVar(randRemoveInds) = NaN;
                
                
                varNames = predictorVarsNames;
                varNames{end+1} = responseVarsNames{ephysInd};
                
                mdl = LinearModel.stepwise(predictorVars(articleInds,:),currRespVar(articleInds),...
                    'linear', 'CategoricalVars', logical(catVarInds), 'VarNames', varNames,...
                    'Criterion', criteria, 'Lower', [varNames{end},' ~ NeuronType'], 'Verbose', 0, 'Upper', complexity,'NSteps',MAXFITSTEPS);
                
%                 bestModelStruct.model{1} = mdl;
                %     mdl.Formula;
                %     bestModelStruct{i}.model{1} = mdl;
                %     mdl
                if strcmp(complexity, 'constant')
                    adjDataMat = ones(size(responseVars,1),1)*nanmean(currRespVar(articleInds));
                else
                    [adjDataMat, actDataMat] = adjustEphysData3(mdl, predictorVarsNames, catVarInds, dataUseInds, responseVars(:,ephysInd));
                end
                currRespVar = responseVars(dataUseInds,ephysInd);
                tempSse = sum((currRespVar(randRemoveInds) - (adjDataMat(randRemoveInds,1))).^2);
                allSses(i,j,k,n) = tempSse;
            end
        end
    end
end

% [neuronCategorical, neuronNames] = grp2idx(NeuronType(dataUseInds));
% numUniqueNeurons = max(neuronCategorical);
% neuronUniqueArtMat = zeros(numUniqueNeurons, numResponseVars);
% for i = 1:numResponseVars
%     currRespVar = responseVars(dataUseInds,i);
%     for j = 1:numUniqueNeurons
%         numUniqueArts = sum(neuronCategorical == j & ~isnan(currRespVar));
%         neuronUniqueArtMat(j,i) = numUniqueArts;
%     end
% end