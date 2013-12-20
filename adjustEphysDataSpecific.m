function [adjDataMat, actDataMat] = adjustEphysDataSpecific(bestModelStruct, predictorVarsNames, catVarInds, dataUseInds, responseVars, specificPredictors)

numEphysVals = length(bestModelStruct);
adjDataMat = [];
actDataMat = [];
for j = 1:numEphysVals
    ephysInd = j;



    neuronOrd = double((bestModelStruct{ephysInd}.model{1}.Variables(:,1)));
    adjVarMat = neuronOrd;

    varNames = get(bestModelStruct{ephysInd}.model{1}.Variables, 'VarNames');

    adjVarVec = [];
    for i = 2:length(varNames)-1
        adjVarVec = [adjVarVec specificPredictors(i)];
    end

    normPredVarMat = [];
    for i = 1:length(varNames)-1
        vec = double(bestModelStruct{ephysInd}.model{1}.Variables(:,i));
        currVarName = varNames{i};
        % if missing predictor value, replace with mode or mean
        if catVarInds(find(strcmp(predictorVarsNames, currVarName)))
            vec(isnan(vec)) = mode(vec(~isnan(vec)));
        else
            vec(isnan(vec)) = nanmean(vec);
        end
        normPredVarMat = vectCat(normPredVarMat, vec);
    end


    obsDataPts = responseVars(dataUseInds,ephysInd);
    fittedDataPts = predict(bestModelStruct{ephysInd}.model{1}, normPredVarMat);


    a = repmat(adjVarVec, length(bestModelStruct{ephysInd}.model{1}.Variables(:,1)), 1);
    adjVarMat = [adjVarMat a];
    tempAdjustment = predict(bestModelStruct{ephysInd}.model{1}, adjVarMat);

    adjVars = (obsDataPts - fittedDataPts) + tempAdjustment;
    actVars = double(bestModelStruct{ephysInd}.model{1}.Variables(:,end));
    
    adjDataMat = vectCat(adjDataMat, adjVars);
    actDataMat = vectCat(actDataMat, actVars);
end