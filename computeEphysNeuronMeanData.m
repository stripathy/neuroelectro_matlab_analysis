function [ephysDataMatNeuronMean, ephysDataMatNeuronStd] = computeEphysNeuronMeanData(ephysDataMat, numUniqueNeurons, neuronUniqueArtMat,neuronCategorical)

% ephysDataMat = adjDataMat;

ephysDataMatNeuronMean = zeros(numUniqueNeurons, size(ephysDataMat,2));
ephysDataMatNeuronStd = zeros(numUniqueNeurons, size(ephysDataMat,2));
for i = 1:numUniqueNeurons
    for j = 1:size(ephysDataMat,2)
        numEphysRepeats = neuronUniqueArtMat(i,j);
        if numEphysRepeats > 0
            ephysDataMatNeuronMean(i,j) = nanmean(ephysDataMat(neuronCategorical==i, j));
            ephysDataMatNeuronStd(i,j) = nanstd(ephysDataMat(neuronCategorical==i, j));
        else 
            ephysDataMatNeuronMean(i,j) = NaN;
            ephysDataMatNeuronStd(i,j) = NaN;
        end
    end
end

