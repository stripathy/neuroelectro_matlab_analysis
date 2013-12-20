function [bestModelStruct, neuronOnlyModelStruct] = metadata_fit_models_criteria(dataInputStruct, fit_complexity_vec, fit_criteria_vec)

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
% fit_complexity = 'purequadratic';
% fit_criteria = 'bic';

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
minIterImprovement = .01;

for i = 1:length(fit_complexity_vec)
    criteria = fit_criteria_vec{i};
    complexity = fit_complexity_vec{i};
%     if i == 3 
%         criteria = 'bic';
%         complexity = fit_complexity_vec;
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
    
%     if sampPct ~= 1
%         articleInds = randsample(articleInds, round(sampPct*length(articleInds)));
%     end
    
    varNames = predictorVarsNames;
    varNames{end+1} = responseVarsNames{ephysInd};
    
    mdl = LinearModel.stepwise(predictorVars(articleInds,:),currRespVar(articleInds),...
            'linear', 'CategoricalVars', logical(catVarInds), 'VarNames', varNames,...
        'Criterion', criteria, 'Lower', [varNames{end},' ~ NeuronType'], 'Verbose', 0, 'Upper', complexity,'NSteps',MAXFITSTEPS);
        
%     mdl = LinearModel.stepwise(predictorVars(articleInds,:),currRespVar(articleInds),...
%         'linear', 'CategoricalVars', logical(catVarInds(1:end)), 'VarNames', varNames,...
%         'Criterion', criteria, 'Lower', [varNames{end},' ~ NeuronType'], 'Verbose', 0, 'Upper', 'quadratic');
    
%     mdl = LinearModel.stepwise(predictorVars(articleInds,2:end),currRespVar(articleInds),...
%         'linear', 'CategoricalVars', logical(catVarInds(2:end)), 'VarNames', varNames(2:end),...
%         'Criterion', criteria, 'Verbose', 0, 'Upper', 'quadratic');
    mdl.Formula;
    bestModelStruct{i}.model{1} = mdl;
%     mdl
    
    neuronOnlyVarNames{1} = varNames{1};
    neuronOnlyVarNames{2} = varNames{end};
    neuronOnlyMdl = LinearModel.fit(predictorVars(articleInds,1),currRespVar(articleInds),...
            'linear', 'CategoricalVars', logical(catVarInds(1)), 'VarNames', neuronOnlyVarNames);
    neuronOnlyModelStruct{i}.model{1} = neuronOnlyMdl;
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