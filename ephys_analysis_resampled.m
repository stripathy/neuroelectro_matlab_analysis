sampPctVec = .2:.1:1;
numSamples = 25;
dataInputStructResamp = [];

cnt = 1;
for o = 1:length(sampPctVec)
    
    sampPct = sampPctVec(o);
    
    for k = 1:numSamples
        responseVarsSampled = responseVars;
        if sampPct ~= 1
            responseVarsSampled = responseVars;
            for i = 1:size(responseVars,1)
                for j = 1:size(responseVars,2)
                    if ~isnan(responseVars(i,j)) && rand > sampPct
                        responseVarsSampled(i,j) = NaN;
                    end
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
        neuronArticleCount = max(neuronUniqueArtMat')';
        
        dataInputStructResamp(cnt).responseVars = responseVarsSampled;
        dataInputStructResamp(cnt).dataUseInds = dataUseInds;
        dataInputStructResamp(cnt).predictorVarsCat = predictorVarsCat;
        dataInputStructResamp(cnt).predictorVarsCont = predictorVarsCont;
        dataInputStructResamp(cnt).predictorVarsNames = predictorVarsNames;
        dataInputStructResamp(cnt).responseVarsNames = responseVarsNames;
        dataInputStructResamp(cnt).catVarInds = catVarInds;
        dataInputStructResamp(cnt).NeuronType = NeuronType;
        dataInputStructResamp(cnt).neuronArticleCount = neuronArticleCount;
        dataInputStructResamp(cnt).neuronUniqueArtMat = neuronUniqueArtMat;
        dataInputStructResamp(cnt).numUniqueNeurons = numUniqueNeurons;
        dataInputStructResamp(cnt).neuronCategorical = neuronCategorical;
        dataInputStructResamp(cnt).sampPct = sampPct;
        cnt = cnt + 1;
    end
end

neuronOnlyModelRsqResamp = zeros(length(dataInputStructResamp), 6);
metadataModelRsqResamp = zeros(length(dataInputStructResamp), 6);
outStruct = [];
parfor k = 1:length(dataInputStructResamp)
    [bestModelStructTemp, neuronOnlyModelStructTemp] = metadata_fit_models_criteria(dataInputStructResamp(k), best_complexity, best_criteria);
    for j = 1:6
        ephysInd = j;
        neuronOnlyModelRsqResamp(k,j) = neuronOnlyModelStructTemp{ephysInd}.model{1}.Rsquared.Adjusted;
        metadataModelRsqResamp(k,j) = bestModelStructTemp{ephysInd}.model{1}.Rsquared.Adjusted;
    end
    outStruct(k).bestModelStruct = bestModelStructTemp;
    outStruct(k).neuronOnlyModelStruct = neuronOnlyModelStructTemp;
end

for k = 1:length(dataInputStructResamp)
    [ephysDistMatResamp, ephysNeuronNamesResamp, ephysWeightsResamp] = compute_neuron_distance_matrix(dataInputStructResamp(k), outStruct(k).bestModelStruct, outStruct(k).neuronOnlyModelStruct, neuronNames);
    outStruct(k).ephysDistMat = ephysDistMatResamp;
    outStruct(k).ephysNeuronNames = ephysNeuronNamesResamp;
    outStruct(k).ephysWeights = ephysWeightsResamp;
end

numSampPts = length(sampPctVec);
% templateDistMat = outStruct(end).ephysDistMat;
% templateNeuronNames = outStruct(end).ephysNeuronNames;
templateDistMat = ephysDistMat;
templateNeuronNames = neuronNames(neuronUseInds(perm));
distMatCorrVec = [];
for k = 1:length(dataInputStructResamp)
    currDistMat = outStruct(k).ephysDistMat;
    currNeuronNames = outStruct(k).ephysNeuronNames;

    allNeuronNames = templateNeuronNames;
    sharedNeuronTypes = {};
    ephysSharedInds = [];
    genesSharedInds = [];
    for i = 1:length(allNeuronNames)
        if  sum(strcmp(allNeuronNames(i), currNeuronNames')) > 0 && sum(strcmp(allNeuronNames(i), templateNeuronNames)) > 0
            genesSharedInds(end+1) = find(strcmp(allNeuronNames(i), currNeuronNames'));
            ephysSharedInds(end+1) = find(strcmp(allNeuronNames(i), templateNeuronNames'));
            sharedNeuronTypes(end+1) = allNeuronNames(i);
        end
    end
    outStruct(k).ephysNeuronInds = genesSharedInds;
    
    ephysSubDistMat = templateDistMat(ephysSharedInds,ephysSharedInds);
    geneSubDistMat = currDistMat(genesSharedInds,genesSharedInds);


    ephysPairwiseDistVals = squareform(ephysSubDistMat)';
    genePairwiseDistVals = squareform(geneSubDistMat)';
    distMatCorrVec(end+1) = corr(ephysPairwiseDistVals, genePairwiseDistVals);
end

distMatCorrMat = reshape(distMatCorrVec, numSamples, numSampPts);
distMatCorrMeanVec = mean(distMatCorrMat);
distMatCorrStdVec = std(distMatCorrMat);


meanEphysSampRsq = zeros(length(sampPctVec), numResponseVars);
stdEphysSampRsq = zeros(length(sampPctVec), numResponseVars);
meanEphysSampRsqNeuron = zeros(length(sampPctVec), numResponseVars);
stdEphysSampRsqNeuron = zeros(length(sampPctVec), numResponseVars);



for i = 1:numResponseVars
    rsqMat = reshape(metadataModelRsqResamp(:,i), numSamples, numSampPts);
    meanEphysSampRsq(:,i) = nanmean(rsqMat);
    stdEphysSampRsq(:,i) = nanstd(rsqMat);
    rsqMat = reshape(neuronOnlyModelRsqResamp(:,i), numSamples, numSampPts);
    meanEphysSampRsqNeuron(:,i) = nanmean(rsqMat);
    stdEphysSampRsqNeuron(:,i) = nanstd(rsqMat);
end