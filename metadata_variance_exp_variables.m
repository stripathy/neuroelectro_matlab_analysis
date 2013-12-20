
% responseVars = [rmp	log10(ir) log10(tau)	amp	hw thresh];
% dataUseInds = log10(Age) > .67 & ~strcmp('cell culture ', PrepType);
% % predictorVarsCat = [NeuronType	Species	Strain	ElectrodeType PrepType];
% predictorVarsCat = [NeuronType	Species	Strain	ElectrodeType PrepType JxnPotential Grandfather];
% % predictorVarsCont = [Age Temp Weight];
% predictorVarsCont = [log10(Age) Temp PubYear];
% 
% predictorVars = [];
% for i = 1:size(predictorVarsCat, 2)
%     predictorVars(:,end+1) = toCategorical(predictorVarsCat(dataUseInds,i));
% end
% 
% for i = 1:size(predictorVarsCont, 2)
%     predictorVars(:,end+1) = predictorVarsCont(dataUseInds, i);
% end
% % predictorVarsNames = {'NeuronType';	'Species';	'Strain';	'ElectrodeType'; 'PrepType'; 'Age'; 'Temp'; 'Weight'};
% predictorVarsNames = {'NeuronType';	'Species';	'Strain';	'ElectrodeType'; 'PrepType'; 'JxnPotential'; 'Grandfather'; 'Age'; 'Temp'; 'PubYear'};
% responseVarsNames = {'rmp'; 'ir'; 'tau'; 'amp'; 'hw'; 'thresh'};
% 
% % catVarInds = [1 1 1 1 1 0 0 0];
% catVarInds = [1 1 1 1 1 1 1 0 0 0];


numResponseVars = size(responseVars, 2);
[neuronCategorical, neuronNames] = grp2idx(NeuronType(dataUseInds));
numUniqueNeurons = max(neuronCategorical);
minNeuronMentions = 1;
minIterImprovement = .01;

for i = [1 2 3 4 5 6]

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
    
    for j = 1:length(predictorVarsNames)
    
        usePredictors = unique([1 j], 'stable');
        varNames = predictorVarsNames(usePredictors);
        varNames{end+1} = responseVarsNames{ephysInd};

        mdl = LinearModel.fit(predictorVars(articleInds,usePredictors),currRespVar(articleInds),...
            'linear', 'CategoricalVars', logical(catVarInds(usePredictors)), 'VarNames', varNames);
        mdl.Formula;
        modelStructAll{i,j}.model{1} = mdl;
    end
end

rsqMat = zeros(size(modelStructAll));

pValMat = zeros(size(modelStructAll));
for i = 1:size(modelStructAll, 1)
    for j = 1:size(modelStructAll, 2)
        rsqMat(i,j) = modelStructAll{i,j}.model{1}.Rsquared.Adjusted;
        tbl = anova(modelStructAll{i,j}.model{1}, 'component');
        pValMat(i,j) = tbl.pValue(end-1);
    end
end
figure; imagesc(rsqMat); 
set(gca, 'XTick', [1:length(predictorVarsNames)], 'XTickLabel', predictorVarsNames);
colormap Hot
xticklabel_rotate;
colorbar;
% set(gca, 'YTick', [1:length(responseVarsNames)], 'YTickLabel', responseVarsNamesPlot);
format_ticks(gca,[], responseVarsNamesPlot)

caxis([0 max(max(rsqMat))]);
hold on;

pValueThresh1 = .05;
pValueThresh2 = .1;
for i = 1:size(modelStructAll, 1)
    for j = 1:size(modelStructAll, 2)
        if (pValMat(i,j) <= pValueThresh1)
            plot(j,i, 'go');
        elseif (pValMat(i,j) <= pValueThresh2)
            plot(j,i, 'g.');
        end
    end
end

rsqDiffMat = rsqMat(:,2:end) - repmat(rsqMat(:,1), 1, size(rsqMat, 2)-1);
figure; imagesc(rsqDiffMat); 
set(gca, 'XTick', [1:length(predictorVarsNames)-1], 'XTickLabel', predictorVarsNames(2:end));
colormap Hot
xticklabel_rotate;
h = colorbar;
ylabel(h, ' Adj. R^2 (relative to neuron-type only model) ')
set(gca, 'YTick', [1:length(responseVarsNames)], 'YTickLabel', responseVarsNames);
format_ticks(gca,[], responseVarsNamesPlot)
caxis([0 max(max(rsqDiffMat))]);
hold on;

for i = 1:size(modelStructAll, 1)
    for j = 2:size(modelStructAll, 2)
        if (pValMat(i,j) <= pValueThresh1)
            plot(j-1,i, 'go');
%         elseif (pValMat(i,j) <= pValueThresh2)
%             plot(j-1,i, 'g.');
        end
    end
end



[neuronCategorical, neuronNames] = grp2idx(NeuronType(dataUseInds));
numUniqueNeurons = max(neuronCategorical);
neuronUniqueArtMat = zeros(numUniqueNeurons, numResponseVars);
for i = 1:numResponseVars
    currRespVar = responseVars(dataUseInds,i);
    for j = 1:numUniqueNeurons
        numUniqueArts = sum(neuronCategorical == j & ~isnan(currRespVar));
        neuronUniqueArtMat(j,i) = numUniqueArts;
    end
end

%% plot specific examples of metadata predictors and ephys props for supp figure

% plot some data that'll go into the main figure
responseVarPlotInds = [2 3 6];
predictorVarPlotInds = [4 7 6];
figure;
for j = 1:length(responseVarPlotInds)
    subplot(1, 3,j);
    
    responseVarInd = responseVarPlotInds(j);
    predictorVarInd = predictorVarPlotInds(j);
    
    tempPredictorVars = double(modelStructAll{responseVarInd,predictorVarInd}.model{1}.Variables(:,1:2));
    currPredictorVars = double(modelStructAll{responseVarInd,predictorVarInd}.model{1}.Variables(:,1:2));
    if catVarInds(predictorVarInd) % check if variable is categorical
        tempPredictorVars(:,2) = mode(tempPredictorVars(:,2));
    else
        tempPredictorVars(:,2) = nanmean(tempPredictorVars(:,2));
    end
    predVals = predict(modelStructAll{responseVarInd,predictorVarInd}.model{1}, tempPredictorVars);
    
    yvals = double(modelStructAll{responseVarInd,predictorVarInd}.model{1}.Variables(:,3)) - predVals + nanmean(responseVars(dataUseInds,responseVarInd));
    if responseVarInd == 2 || responseVarInd == 5 || responseVarInd == 3
        yvals = 10.^yvals;
    end
    if catVarInds(predictorVarInd) % check if variable is categorical
        metadataVals = predictorVarsCat(dataUseInds, predictorVarInd);
        numMetadataVals = length(unique(metadataVals));
        metadataXVals = unique(metadataVals, 'stable');
        
        uniquePredictorVars = unique(currPredictorVars(:, 2), 'stable');
        numUniquePredictorVars = numel(uniquePredictorVars);
        predMeans = [];
        predStds = [];
        for i = 1:numUniquePredictorVars
            currMetadataInds = currPredictorVars(:, 2) == i;
            predMeans(i) = nanmean(yvals(currMetadataInds));
            predStds(i) = nanstd(yvals(currMetadataInds));
        end
        
        
        jitterInds = randn(length(yvals), 1)*.05;
        
        if responseVarInd == 2 || responseVarInd == 3

            semilogy(currPredictorVars(:, 2)+jitterInds, yvals, 'k.')
            hold on;
            set(gca, 'XTick', uniquePredictorVars, 'XTickLabel', metadataXVals);

            errorbar(uniquePredictorVars+.25, predMeans, predStds,'k.');
            semilogy(uniquePredictorVars+.25, predMeans, 'k*');
        else
            plot(currPredictorVars(:, 2)+jitterInds, yvals, 'k.')
            hold on;
            set(gca, 'XTick', uniquePredictorVars, 'XTickLabel', metadataXVals);

            errorbar(uniquePredictorVars+.25, predMeans, predStds,'k.');
            plot(uniquePredictorVars+.25, predMeans, 'k*');
        end
            box off;
            ylabel([ephysVarsNames{responseVarInd} ' (' ephysVarsUnits{responseVarInd} ') ']);
    else
%         subplot(2,2,j);
        jitterInds = randn(length(yvals), 1)*.05;
        
        if predictorVarInd == 7
            if responseVarInd == 2 
                loglog(10.^(currPredictorVars(:, 2)+jitterInds), yvals, 'k.')
            else
                semilogx(10.^(currPredictorVars(:, 2)+jitterInds), yvals, 'k.')
            end
        else
            plot(currPredictorVars(:, 2)+jitterInds, yvals, 'k.')
        end
        hold on;
        box off;
        ylabel([ephysVarsNames{responseVarInd} ' (' ephysVarsUnits{responseVarInd} ') ']);
        if predictorVarInd == 8
            xlabel(' Recording temperature (°C) ');
            axis([20 38 0 max(yvals)*1.1]);
        elseif predictorVarInd == 7
            xlabel(' Animal age (days) ');
            set(gca, 'XTick', [10 100])
            axis([5 500 0 max(yvals)*1.1]);
        end
    end
end

% plot some supp figures

responseVarPlotInds = [2 6 4 5 3 6];
predictorVarPlotInds = [4 4 5 8 8 6];
figure;
for j = 1:length(responseVarPlotInds)
    subplot(3,2,j);
    
    responseVarInd = responseVarPlotInds(j);
    predictorVarInd = predictorVarPlotInds(j);
    
    tempPredictorVars = double(modelStructAll{responseVarInd,predictorVarInd}.model{1}.Variables(:,1:2));
    currPredictorVars = double(modelStructAll{responseVarInd,predictorVarInd}.model{1}.Variables(:,1:2));
    if catVarInds(predictorVarInd) % check if variable is categorical
        tempPredictorVars(:,2) = mode(tempPredictorVars(:,2));
    else
        tempPredictorVars(:,2) = nanmean(tempPredictorVars(:,2));
    end
    predVals = predict(modelStructAll{responseVarInd,predictorVarInd}.model{1}, tempPredictorVars);
    
    yvals = double(modelStructAll{responseVarInd,predictorVarInd}.model{1}.Variables(:,3)) - predVals + nanmean(responseVars(dataUseInds,responseVarInd));
    if responseVarInd == 2 || responseVarInd == 5 || responseVarInd == 3
        yvals = 10.^yvals;
    end
    if catVarInds(predictorVarInd) % check if variable is categorical
        metadataVals = predictorVarsCat(dataUseInds, predictorVarInd);
        numMetadataVals = length(unique(metadataVals));
        metadataXVals = unique(metadataVals, 'stable');
        
        uniquePredictorVars = unique(currPredictorVars(:, 2), 'stable');
        numUniquePredictorVars = numel(uniquePredictorVars);
        predMeans = [];
        predStds = [];
        for i = 1:numUniquePredictorVars
            currMetadataInds = currPredictorVars(:, 2) == i;
            predMeans(i) = nanmean(yvals(currMetadataInds));
            predStds(i) = nanstd(yvals(currMetadataInds));
        end
        
        
        jitterInds = randn(length(yvals), 1)*.05;
        
        if responseVarInd == 2

            semilogy(currPredictorVars(:, 2)+jitterInds, yvals, 'k.')
            hold on;
            set(gca, 'XTick', uniquePredictorVars, 'XTickLabel', metadataXVals);

            errorbar(uniquePredictorVars+.25, predMeans, predStds,'k.');
            semilogy(uniquePredictorVars+.25, predMeans, 'k*');
        else
            plot(currPredictorVars(:, 2)+jitterInds, yvals, 'k.')
            hold on;
            set(gca, 'XTick', uniquePredictorVars, 'XTickLabel', metadataXVals);

            errorbar(uniquePredictorVars+.25, predMeans, predStds,'k.');
            plot(uniquePredictorVars+.25, predMeans, 'k*');
        end
            box off;
            ylabel([ephysVarsNames{responseVarInd} ' (' ephysVarsUnits{responseVarInd} ') ']);
    else
%         subplot(2,2,j);
        jitterInds = randn(length(yvals), 1)*.05;
        
        if predictorVarInd == 7
            semilogx(10.^(currPredictorVars(:, 2)+jitterInds), yvals, 'k.')
        else
            plot(currPredictorVars(:, 2)+jitterInds, yvals, 'k.')
        end
        hold on;
        box off;
        ylabel([ephysVarsNames{responseVarInd} ' (' ephysVarsUnits{responseVarInd} ') ']);
        if predictorVarInd == 8
            xlabel(' Recording temperature (°C) ');
            axis([20 38 0 max(yvals)*1.1]);
        elseif predictorVarInd == 7
            xlabel(' Animal age (days) ');
            set(gca, 'XTick', [10 100])
            axis([5 500 0 max(yvals)*1.1]);
        end
    end
end

