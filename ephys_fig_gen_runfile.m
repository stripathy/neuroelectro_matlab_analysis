%% load ephys data from spreadsheet
ephys_metadata_import_fxn

%% plot basic scatter points of ephys data
ephys_desc_fig_gen

%% plot histograms showing spread in expt metadata
metadata_histo_examples_fig_gen

%% generate regression models between metadata and ephys
% set up data structure

% perform basic data filtering - 
% filter out neuron types not in brain (these weren't manually validated so closely)
% filter out time constant values diff from 2sigma from median

dataUseInds = log10(Age) > .67 & ~strcmp('cell culture ', PrepType) ;
speciesUseInds = strcmp('Rats ', Species) | strcmp('Mice ', Species) | strcmp('Guinea Pigs ', Species);
dataUseInds = dataUseInds & ~strcmp('Neocortex uncharacterized cell ', NeuronType);
dataUseInds = dataUseInds & ~strcmp('Cochlea hair cell inner ', NeuronType);
dataUseInds = dataUseInds & ~strcmp('Dorsal root ganglion cell ', NeuronType);
dataUseInds = dataUseInds & ~strcmp('DRG temperature cell ', NeuronType);
dataUseInds = dataUseInds & ~strcmp('Spinal cord intermediate horn motor neuron sympathetic ', NeuronType);
dataUseInds = dataUseInds & ~strcmp('Spinal cord ventral horn interneuron IA ', NeuronType);
dataUseInds = dataUseInds & ~strcmp('Spinal cord ventral horn interneuron II ', NeuronType);
dataUseInds = dataUseInds & ~strcmp('Spinal cord ventral horn interneuron V2 ', NeuronType);
dataUseInds = dataUseInds & ~strcmp('Spinal cord ventral horn motor neuron alpha ', NeuronType);

dataUseInds = dataUseInds & speciesUseInds;

responseVars = [rmp	log10(ir) tau	amp	hw thresh];
thresholdsd = 2;
badVarIndsAll = [];
inds1 = find(dataUseInds);
for i = [3]
currResponseVars = responseVars(dataUseInds,i);
tempVarMed = nanmedian(currResponseVars);
tempVarStd = nanstd(currResponseVars);
badVarInds = find(currResponseVars < (tempVarMed - thresholdsd*tempVarStd) |...
    currResponseVars > (tempVarMed + thresholdsd*tempVarStd));
badVarIndsAll = vectCat(badVarIndsAll, badVarInds);
responseVars(inds1(badVarInds),i) = NaN;
%     responseVars(,i) = NaN;
end
responseVars(:,3) = log10(responseVars(:,3));
responseVars(:,5) = log10(responseVars(:,5));

% responseVars = [rmp	log10(ir) log10(tau)	amp	log10(hw) thresh];


%% set up basic data structure to use for metadata prediction modeling
predictorVarsCat = [NeuronType	Species	Strain	ElectrodeType PrepType JxnPotential];
predictorVarsCont = [log10(Age) Temp];

predictorVarsNames = {'NeuronType';	'Species';	'Strain';	'ElectrodeType'; 'PrepType'; 'JxnPotential'; 'Age'; 'Temp'};
responseVarsNames = {'rmp'; 'ir'; 'tau'; 'amp'; 'hw'; 'thresh'};

predictorVars = [];
for i = 1:size(predictorVarsCat, 2)
    predictorVars(:,end+1) = toCategorical(predictorVarsCat(dataUseInds,i));
end

for i = 1:size(predictorVarsCont, 2)
    predictorVars(:,end+1) = predictorVarsCont(dataUseInds, i);
end

catVarInds = [1 1 1 1 1 1 0 0];

responseVarsSampled = responseVars;
    
clear dataInputStruct
dataInputStruct.responseVars = responseVarsSampled;
dataInputStruct.dataUseInds = dataUseInds;
dataInputStruct.predictorVarsCat = predictorVarsCat;
dataInputStruct.predictorVarsCont = predictorVarsCont;
dataInputStruct.predictorVarsNames = predictorVarsNames;
dataInputStruct.responseVarsNames = responseVarsNames;
dataInputStruct.catVarInds = catVarInds;
dataInputStruct.NeuronType = NeuronType;

%% fit best models using all data based on cross val results from metadata_crossval_run

% these best complexity and criteria measures come from cross-validation
% (analyses below)
best_complexity = {'linear','linear','linear','linear','linear','purequadratic'};
best_criteria = {'bic', 'bic', 'bic', 'bic', 'bic', 'bic'};

% runs the analysis
[bestModelStruct, neuronOnlyModelStruct] = metadata_fit_models_criteria(dataInputStruct, best_complexity, best_criteria);



%% plot data example showing importance of metadata

% plot just data from ca1 and msns

ephysInd = 2;
minNeuronMentions = 1;
currRespVar = responseVars(dataUseInds,ephysInd);

numResponseVars = size(responseVars, 2);
[neuronCategorical, neuronNames] = grp2idx(NeuronType(dataUseInds));

m = find(strcmp('Hippocampus CA1 pyramidal cell ', neuronNames));
% m = find(strcmp('Neocortex basket cell', neuronNames));

articleInds = find(neuronCategorical == m & ~isnan(currRespVar));

patchInds = find(predictorVars(articleInds,4) == 1);
sharpInds = find(predictorVars(articleInds,4) == 2);

figure; loglog(10.^predictorVars(articleInds(patchInds),7), 10.^currRespVar(articleInds(patchInds)), 'bo')

hold on;
loglog(10.^predictorVars(articleInds(sharpInds),7), 10.^currRespVar(articleInds(sharpInds)), 'b*')


m = find(strcmp('Neostriatum medium spiny neuron ', neuronNames));
articleInds = find(neuronCategorical == m & ~isnan(currRespVar));
patchInds = find(predictorVars(articleInds,4) == 1);
sharpInds = find(predictorVars(articleInds,4) == 2);

loglog(10.^predictorVars(articleInds(patchInds),7), 10.^currRespVar(articleInds(patchInds)), 'go')
loglog(10.^predictorVars(articleInds(sharpInds),7), 10.^currRespVar(articleInds(sharpInds)), 'g*')
xlabel('Animal Age (days)');
ylabel('Input Resistance (MOhms)');
axis([10 200 10 1000])
box off;
%% plot set of bar graphs showing adj rsq for models with and without metadata adjusting

neuronOnlyModelRsq = [];
metadataModelRsq = [];
pValueMetadataTerms = [];

for j = 1:6
    ephysInd = j;
    neuronOnlyModelRsq(j) = neuronOnlyModelStruct{ephysInd}.model{1}.Rsquared.Adjusted;
    metadataModelRsq(j) = bestModelStruct{ephysInd}.model{1}.Rsquared.Adjusted;
    
    neuronOnlyModelTermsCount = length(neuronOnlyModelStruct{ephysInd}.model{1}.CoefficientNames);
    metadataModelTermsCount = length(bestModelStruct{ephysInd}.model{1}.CoefficientNames);
    termsVec = zeros(1, metadataModelTermsCount);
    termsVec(neuronOnlyModelTermsCount+1:metadataModelTermsCount) = 1;
    %pValueMetadataTerms(j) = coefTest(bestModelStruct{ephysInd}.model{1},termsVec);
end

barColors = ['k', 'r'];
figure; hold on;
barInd = 2;
for j = 1:6
    bar(barInd, neuronOnlyModelRsq(j), barColors(1));
    barInd = barInd + 1;
    bar(barInd, metadataModelRsq(j), barColors(2));
    barInd = barInd + 2;
end
ylabel('Model adj. R^2')
% set(gca, 'XTick', 2:3:19, 'XTickLabel', responseVarsNames)
format_ticks(gca,responseVarsNamesPlot, [], 3:3:20, [], [])
set(gca, 'XTick', [])

%% plot influence of each metadata predictor on each ephys term
metadata_variance_exp_variables

%% adjust ephys data - but to a specific criterion set
predictorVarsNames = {'NeuronType';	'Species';	'Strain';	'ElectrodeType'; 'PrepType'; 'JxnPotential'; 'Age'; 'Temp'};

% these are what ephys metadata conditions used in the UrbanLab data
predictorVarsSpecificNames = {NeuronType{1}, 'Mice ', 'C57BL ', 'Patch-clamp ','in vitro ', 'Not reported ', log10(20), 37};

predictorVarsSpecific = [];

for i = 1:length(predictorVarsSpecificNames)
    if catVarInds(find(strcmp(predictorVarsNames, predictorVarsNames{i}))) % is categorical
        currPredVars = unique(predictorVarsCat(dataUseInds,i), 'stable');
        idx = find(strcmp(currPredVars, predictorVarsSpecificNames{i}));
        predictorVarsSpecific(i) = idx;
    else
        predictorVarsSpecific(i) = predictorVarsSpecificNames{i};
    end
end
%     index = strfind(currPredVars, predictorVarsSpecific{i})

[adjDataMatSpecific, actDataMat] = adjustEphysDataSpecific(bestModelStruct, predictorVarsNames, catVarInds, dataUseInds, responseVars,predictorVarsSpecific);

% diff between Specific and not is that data has been un-log transformed
actDataMat(:,2) = 10.^actDataMat(:,2);
actDataMat(:,3) = 10.^actDataMat(:,3);
actDataMat(:,5) = 10.^actDataMat(:,5);
adjDataMatSpecific(:,2) = 10.^adjDataMatSpecific(:,2);
adjDataMatSpecific(:,3) = 10.^adjDataMatSpecific(:,3);
adjDataMatSpecific(:,5) = 10.^adjDataMatSpecific(:,5);

[ephysMeanDataSpecific, ephysMeanDataStdSpecific] = computeEphysNeuronMeanData(adjDataMatSpecific, numUniqueNeurons, neuronUniqueArtMat,neuronCategorical);
[ephysMeanDataActual, ephysMeanDataStdActual] = computeEphysNeuronMeanData(actDataMat, numUniqueNeurons, neuronUniqueArtMat,neuronCategorical);

outputNeuronNames = NeuronType(dataUseInds);
%% plot ca1 and msn ir's adjusted for metadata


% adjusts ephys data for entire dataset to the mean and mode across
% metadata attributes
% make sure neuron type info in spreadsheet is sorted!!
[adjDataMat, actDataMat] = adjustEphysData(bestModelStruct, predictorVarsNames, catVarInds, dataUseInds, responseVars);

ephysInd = 2;
ephysValsAct = actDataMat(:,ephysInd);
ephysValsAdj = adjDataMat(:,ephysInd);

ca1Inds = find(neuronCategorical == find(strcmp('Hippocampus CA1 pyramidal cell ', neuronNames)));
msnInds = find(neuronCategorical == find(strcmp('Neostriatum medium spiny neuron ', neuronNames)));

figure;hold on;  
jitterInds1 = repmat(1, length(ca1Inds), 1) + randn(length(ca1Inds), 1)*.1;
jitterInds2 = repmat(1, length(ca1Inds), 1) + jitterInds1;
plot(jitterInds1, 10.^ephysValsAct(ca1Inds), 'ko');
plot(jitterInds2, 10.^ephysValsAdj(ca1Inds), 'ro');
% errorbar(1-.4, 10.^nanmean(ephysValsAct(ca1Inds)), nanstd(10.^ephysValsAct(ca1Inds)),'k.');
plot(1-.4, 10.^nanmean(ephysValsAct(ca1Inds)), 'k*');
% errorbar(2+.4, 10.^nanmean(ephysValsAdj(ca1Inds)), nanstd(10.^ephysValsAdj(ca1Inds)),'r.');
plot(2+.4, 10.^nanmean(ephysValsAdj(ca1Inds)), 'r*');

jitterInds1 = repmat(4, length(msnInds), 1) + randn(length(msnInds), 1)*.1;
jitterInds2 = repmat(1, length(msnInds), 1) + jitterInds1;
plot(jitterInds1, 10.^ephysValsAct(msnInds), 'ko');
plot(jitterInds2, 10.^ephysValsAdj(msnInds), 'ro');
% errorbar(4-.4, 10.^nanmean(ephysValsAct(msnInds)), nanstd(10.^ephysValsAct(msnInds)),'k.');
plot(4-.4, 10.^nanmean(ephysValsAct(msnInds)), 'k*');
% errorbar(5+.4, 10.^nanmean(ephysValsAdj(msnInds)), nanstd(10.^ephysValsAdj(msnInds)),'r.');
plot(5+.4, 10.^nanmean(ephysValsAdj(msnInds)), 'r*');

errorbar([1-.4 4-.4], [10.^nanmean(ephysValsAct(ca1Inds)) 10.^nanmean(ephysValsAct(msnInds))], ...
    [nanstd(10.^ephysValsAct(ca1Inds)) nanstd(10.^ephysValsAct(msnInds))], 'k.');

errorbar([2+.4 5+.4], [10.^nanmean(ephysValsAdj(ca1Inds)) 10.^nanmean(ephysValsAdj(msnInds))], ...
    [nanstd(10.^ephysValsAdj(ca1Inds)) nanstd(10.^ephysValsAdj(msnInds))], 'r.');


axis([0 6 0 350])
ylabel('Input resistance (MOhms)')
set(gca, 'XTick', [1.5 4.5], 'XTickLabel', {' CA1, Pyr '; ' Str, MSN '})

%% plot ephys properties correlation figure

ephysDataMat = adjDataMat;

ephysDataMatNeuronMean = zeros(numUniqueNeurons, size(ephysDataMat,2));
for i = 1:numUniqueNeurons
    for j = 1:size(ephysDataMat,2)
        numEphysRepeats = neuronUniqueArtMat(i,j);
        if numEphysRepeats > 0
            ephysDataMatNeuronMean(i,j) = nanmean(ephysDataMat(neuronCategorical==i, j));
            ephysDataMatNeuronStd(i,j) = nanstd(ephysDataMat(neuronCategorical==i, j));
        else 
            ephysDataMatNeuronMean(i,j) = NaN;
        end
    end
end

ephysData = ephysDataMatNeuronMean;

neuronArticleCount = max(neuronUniqueArtMat')';
neuronMinCount = 2;
% calculate pairwise correlations for ephys values

corrMat = zeros(size(ephysData,2), size(ephysData,2));
pValMat = ones(size(ephysData,2), size(ephysData,2));
for i = 1:size(ephysData,2)
    for j = 1:size(ephysData,2)
        validInds = ~isnan(ephysData(:,i)) & ~isnan(ephysData(:,j));
        if sum(validInds>0) >= 10
            [r_spearman, p_spearman] = corr(ephysData(validInds,i), ephysData(validInds,j), 'type', 'Spearman');
            [r_pearson, p] = corr(ephysData(validInds,i), ephysData(validInds,j), 'type', 'Pearson');
%             [rMat, pMat] = corrcoef(ephysData(validInds,i), ephysData(validInds,j));
%             corrMat(i,j) = rMat(1,2);
%             pValMat(i,j) = pMat(1,2);
             corrMatSpearman(i,j) = r_spearman;
               corrMatPearson(i,j) = r_pearson;
               
               pValMat(i,j) = p_spearman;
               if i == j
                   pValMat(i,j) = 1;
               end
        end
    end
end
corrMat = corrMatSpearman;

% numEphysProps = size(ephysList,1);
% corrMat = corrMat + corrMat';
% pValMat = pValMat + pValMat';
[row, col] = find(pValMat <= .05);
ephysInds = 1:6;

ephysCorrMat = corrMat(ephysInds, ephysInds);
figure; 
% subplot(2,2,3);
% colormap(redblue);
imagesc(log10(pValMat));  

hC = colorbar;
colormap Hot;
caxis([-3 0]);
L = [.001 .005 0.01 0.02 0.05 0.1 0.2 0.5 1];
% Choose appropriate
% or somehow auto generate colorbar labels
l = log10(L); % Tick mark positions
set(hC,'Ytick',l,'YTicklabel',L);


%caxis([-1 1]);
% 
% h = colorbar;
ylabel(hC, ' Correlation p-value ')
hold on;
plot(col, row, 'wo')
% set(gca,'XDir','normal')
format_ticks(gca,responseVarsNamesPlot, responseVarsNamesPlot, [], [], 90, []);
set(gca,'YDir','reverse')

figure; 
% subplot(2,2,3);
colormap(redblue);
imagesc(ephysCorrMat);  caxis([-1 1]);
% 
h = colorbar;
ylabel(h, ' Correlation (Spearman) ')
hold on;
plot(col, row, 'wo')
% set(gca,'XDir','normal')
format_ticks(gca,responseVarsNamesPlot, responseVarsNamesPlot, [], [], 90, []);
set(gca,'YDir','reverse')
% set(gca, 'XTick', [])
% set(gca, 'YTick', 1:length(ephysInds), 'YTickLabel', responseVarsNames(ephysInds))
% set(gca, 'XTick', 1:length(ephysInds), 'XTickLabel', responseVarsNames(ephysInds))
% xticklabel_rotate(1:length(ephysInds),90,responseVarsNames(ephysInds),'interpreter','none')

% plot a couple property pairs for giggles
figure;
subplot(2,2,2);
prop1 = 6;
prop2 = 5;
validInds = ~isnan(ephysData(:,prop1)) & ~isnan(ephysData(:,prop2));

% p = polyfit(ephysData(validInds,prop1),10.^ephysData(validInds,prop2),1);
% fit_points = polyval(p,ephysData(validInds,prop1));
% plot([min(ephysData(validInds,prop1)) max(ephysData(validInds,prop1))], [min(fit_points) max(fit_points)], '-k');
% 
% hold on;

plot(ephysData(validInds,prop1), 10.^ephysData(validInds,prop2) , 'o');
xlabel([ephysVarsNames{prop1} ' (' ephysVarsUnits{prop1} ') ']);
ylabel([ephysVarsNames{prop2} ' (' ephysVarsUnits{prop2} ') ']);



box off;
axis([-55 -30 0 2.5]);

% figure;
% prop1 = 2;
% prop2 = 1;
% validInds = ~isnan(ephysData(:,prop1)) & ~isnan(ephysData(:,prop2));
% semilogx(10.^ephysData(validInds,prop1), ephysData(validInds,prop2) , 'o');
% xlabel([ephysVarsNames{prop1} ' (' ephysVarsUnits{prop1} ') ']);
% ylabel([ephysVarsNames{prop2} ' (' ephysVarsUnits{prop2} ') ']);
% box off;


% figure;
subplot(2,2,1);
prop1 = 2;
prop2 = 3;
validInds = ~isnan(ephysData(:,prop1)) & ~isnan(ephysData(:,prop2));
% p = polyfit(10.^ephysData(validInds,prop1),10.^ephysData(validInds,prop2),1);
% fit_points = polyval(p,10.^ephysData(validInds,prop1));
% semilogx([min(10.^ephysData(validInds,prop1)) max(10.^ephysData(validInds,prop1))], [min(fit_points) max(fit_points)], '-k');
% hold on;
semilogx(10.^ephysData(validInds,prop1), 10.^ephysData(validInds,prop2) , 'o');
xlabel([ephysVarsNames{prop1} ' (' ephysVarsUnits{prop1} ') ']);
ylabel([ephysVarsNames{prop2} ' (' ephysVarsUnits{prop2} ') ']);
box off;
axis([30 1000 0 60]);

% figure;
% prop1 = 3;
% prop2 = 5;
% validInds = ~isnan(ephysData(:,prop1)) & ~isnan(ephysData(:,prop2));
% plot(10.^ephysData(validInds,prop1), 10.^ephysData(validInds,prop2) , 'o');
% xlabel([ephysVarsNames{prop1} ' (' ephysVarsUnits{prop1} ') ']);
% ylabel([ephysVarsNames{prop2} ' (' ephysVarsUnits{prop2} ') ']);
% box off;

% plot pca decomp of matrix
[coeff, latent, pcaExplained] = pcacov(corrMatPearson);



% figure;

% bar(pcaExplained,'k')
% ylabel('Variance explained (%)');
% xlabel('Principal components');
% box off;
% axis([0 7 0 max(pcaExplained)*1.1])

% 
% ca1CorrMat = zeros(size(ephysData,2), size(ephysData,2));
% ca1PValMat = ones(size(ephysData,2), size(ephysData,2));
% ca1Inds = (neuronCategorical==26);
% for i = 1:size(ephysData,2)
%     for j = 1:size(ephysData,2)
%         validInds = ~isnan(ephysDataMat(:,i)) & ~isnan(ephysDataMat(:,j)) & ca1Inds;
%         if sum(validInds>0) >= 3
%             [r, p] = corr(ephysDataMat(validInds,i), ephysDataMat(validInds,j), 'type', 'Pearson');
% %             [rMat, pMat] = corrcoef(ephysData(validInds,i), ephysData(validInds,j));
% %             corrMat(i,j) = rMat(1,2);
% %             pValMat(i,j) = pMat(1,2);
%                ca1CorrMat(i,j) = r;
%                ca1PValMat(i,j) = p;
%                if i == j
%                    ca1PValMat(i,j) = 1;
%                end
%         end
%     end
% end
% [~, ~, ca1PCAExplained] = pcacov(ca1CorrMat);
% 
% [row, col] = find(ca1PValMat <= .05);
% ephysInds = 1:6;
% 
% ephysCorrMat = ca1CorrMat(ephysInds, ephysInds);
% figure; 
% colormap(redblue);
% imagesc(ephysCorrMat); colorbar; caxis([-1 1]);
% hold on;
% plot(col, row, 'wo')
% set(gca, 'YTick', 1:length(ephysInds), 'YTickLabel', responseVarsNames(ephysInds))
% set(gca, 'XTick', 1:length(ephysInds), 'XTickLabel', responseVarsNames(ephysInds))
% xticklabel_rotate(1:length(ephysInds),90,responseVarsNames(ephysInds),'interpreter','none')

%% compute pcas of ephys data (using probabilistic pca)
[neuronCategorical, neuronNames] = grp2idx(NeuronType(dataUseInds));
% compute ephysWeights
ephysWeights = zeros(length(ephysInds), 1);
for j = 1:6
    ephysWeights(j) = metadataModelRsq(j)./max(metadataModelRsq);
end

% compute ppcas of ephys space (note that the algorithm is stochastic, so
% diff results on runs are expected)

% only use neurons where no more than 2 (of the 6 ephys params missing and
% at least 3 articles)
maxNonZeros = 2; % how many ephys params missing?
neuronMinCount = 3; % how many articles defining neuron type?


%runs ppca algorithm
ephysUseInds = 1:6;
neuronUseInds = find(sum(isnan(ephysDataMatNeuronMean(:,ephysUseInds)),2) <= maxNonZeros & neuronArticleCount >= neuronMinCount);
useMat = ephysDataMatNeuronMean(neuronUseInds, ephysUseInds);
for i = 1:size(useMat,2)
    useMat(:,i) = (useMat(:,i) - nanmean(useMat(:,i)))./nanstd(useMat(:,i)) .* ephysWeights(i);
end
[coeff,score,pcvar, mu, v, S] = ppca(useMat, length(ephysUseInds) - 1);
recon = S.Recon;
keepData = recon;

figure;
hold on;
for i = 1:length(score)
    if neuronArticleCount(neuronUseInds(i)) >= 1
        plot(-score(i,1), score(i,2), '.')
        text(-score(i,1), score(i,2)+.001, neuronNames{neuronUseInds(i)});
    end
end
box off;
% axis tight;
ylabel('Principal component 2');
xlabel('Principal component 1');
axis([-2.1 4 -2.5 2]);

% assign neuron 2d prin comp data

figure;
hold on;
for i = 1:length(score)
    if neuronArticleCount(neuronUseInds(i)) >= 1
        plot(find(sort(score(:,1))==score(i,1)), find(sort(score(:,2))==score(i,2)), '.')
        text(find(sort(score(:,1))==score(i,1)), find(sort(score(:,2))==score(i,2))+1, neuronNames{neuronUseInds(i)});
    end
end
box off;
% axis tight;
ylabel('Principal component 2');
xlabel('Principal component 1');

clear nNames xInds yInds
for i = 1:length(score)
    nNames{i,1} = neuronNames{neuronUseInds(i)};
    xInds(i,1) = find(sort(score(:,1))==score(i,1))./length(score(:,1));
    yInds(i,1) = find(sort(score(:,2))==score(i,2))./length(score(:,2));
end

pcaExplained = pcvar./(sum(pcvar))*100;
% figure; 
% figure;
subplot(2,2,4);
bar(pcaExplained,'k')
ylabel('Variance explained (%)');
xlabel('Principal components');
box off;
axis([0 6 0 max(pcaExplained)*1.1])

figure;
bar(-coeff(:,1), 'k');
set(gca,'YDir','normal')
format_ticks(gca,responseVarsNamesPlot, [], [], [], 90, []);
box off;
% ylabel('Variance explained (%)');
% set(gca,'YDir','reverse')







%% plot clustergram of ephys data based on ppcas
metric = 'euclidean';
linkageMethod = 'ward';

ephysUseInds = 1:6;
cg = clustergram(keepData, 'RowLabels', neuronNames(neuronUseInds), 'ColumnLabels', responseVarsNames(ephysUseInds), ...
    'Cluster', 3, 'Standardize', 'none', 'OptimalLeafOrder', 'false', 'Colormap', 'redblue', ...
    'ColumnPDist', metric, 'RowPDist', metric, 'Linkage', linkageMethod);



% compute surrogate samples by drawing distance matrix from ephys corr
% matrix
% cholMat = chol(corrMatPearson);
% numRandSamples = 1000;
% distPVal = .05;
% distPVal2 = .1;
% obs = cholMat*randn(6, numRandSamples);
% obs = obs';
% obs = obs .* repmat(std(recon), size(obs, 1), 1);
% distVecNull = pdist(obs, metric);
% sortedNullDists = sort(distVecNull, 'ascend');
% pValDist1 = sortedNullDists(distPVal*length(sortedNullDists));
% pValDist2 = sortedNullDists(distPVal2*length(sortedNullDists));


% compute distance matrix reflecting clustergram

[neuronCategorical, neuronNames] = grp2idx(NeuronType(dataUseInds));

Z = linkage(keepData, linkageMethod, metric);
figure;
[H, t, perm] = dendrogram(Z,length(score), 'labels', neuronNames(neuronUseInds),   'Orientation', 'left');
xticklabel_rotate
perm = perm(end:-1:1);
% hold on;
% plot([pValDist1 pValDist1],[0 max(perm)+1], 'g--')


ephysUseInds = [2 3 5 6 1 4];
permRevOrder = perm(end:-1:1);
figure; imagesc(keepData(permRevOrder, ephysUseInds));
colormap 'redblue'
colorbar;
set(gca, 'YTick', 1:length(keepData), 'YTickLabel', neuronNames(neuronUseInds(permRevOrder)));
set(gca, 'XTick', 1:length(responseVarsNames(ephysUseInds)), 'XTickLabel', responseVarsNames(ephysUseInds));
% xticklabel_rotate
colorbar;

distMat = squareform(pdist(keepData, metric));
figure;
imagesc(distMat(perm, perm));
colormap('redblue')

h = colorbar;
ylabel(h, ' Electrophysiological distance (Euclidean) ')
set(gca, 'YTick', 1:length(keepData), 'YTickLabel', neuronNames(neuronUseInds(perm)));
set(gca, 'XTick', 1:length(keepData), 'XTickLabel', neuronNames(neuronUseInds(perm)));


xticklabel_rotate

ephysDistMat = distMat(perm, perm);
ephysNeuronNames = neuronNames(neuronUseInds(perm));
% 
% figure;
% [nullDistCDF, nullDistX] = ecdf(sortedNullDists);
% [empDistCDF, empDistX] = ecdf(squareform(distMat));
% % 
% plot(empDistX, empDistCDF, 'k', nullDistX, nullDistCDF, 'r')
% hold on;
% plot([pValDist1 pValDist1], [0 1], 'g--');
% plot([0 max(nullDistX)], [distPVal distPVal], 'g--');
% ylabel('Empirical CDF (Probability)');
% xlabel('Electrophysiological pairwise distance (Euclidean)');
% box off;
% axis tight;
% 
% % 
% binEdges = linspace(min(squareform(distMat)), max(squareform(distMat)), 25);
% [nullDistCDF, nullDistX] = histc(sortedNullDists,binEdges);
% [empDistCDF, empDistX] = histc(squareform(distMat),binEdges);
% 
% plot(binEdges, empDistCDF./sum(empDistCDF), 'k.', binEdges, nullDistCDF./sum(nullDistCDF), 'r.')
% 
% plot(binEdges, empDistCDF./sum(empDistCDF), 'k')
% hold on;
% plot([pValDist1 pValDist1], [0 max(empDistCDF./sum(empDistCDF))], 'g--')


%% load gene expression data (go get a coffee cause it takes a while)
gene_expression_import
load('data/go_data_6_19_13')
load('data/non_exp_allen_genes');

go_gene_map_full(0) = keys(gene_ise_index_map);
go_gene_map_full(-1) = nonExpGenes;

goStruct = GO.terms;
go_id_term_map = containers.Map('KeyType', 'double', 'ValueType','any');
for i = 1:length(goStruct)
    go_id_term_map(goStruct(i).id) = goStruct(i).name;
end
go_id_term_map(0) = 'all genes';
go_id_term_map(-1) = 'nonexpressed genes';

% maps neuron names to neuron regions
validNeuronInds = [];
neuronNames = {};
numNeurons = length(neuronStruct);
for i = 1:numNeurons
    if isfield(neuronStruct{i}, 'regionListIdsPerm')
        validNeuronInds(end+1) = i;
    end
    neuronNames{i} = neuronStruct{i}.name;
    if isfield(neuronStruct{i}, 'regionListNames')
        regionNames{i} = neuronStruct{i}.regionListNames{1};
    end
end

%% produce gene expression, ephys paired distance matrices for ion channel genes

% assign go id, 0 = all genes, 5216 = GO id for ion channel activity
go_id = 5216;

mantel_resamples = 1000;
gene_neuron_correlation_test

%% run code for all GO classes

% you should turn parallel matlab first
% matlabpool local
gene_ephys_correlate_gos_all

%% analyze GO class stuff

[s, sortInds] = sort(corrVals, 'descend');

goIds = checkGos(sortInds);
bestGoTerms = {};
for i = 1:length(goIds)
    bestGoTerms{i} = go_id_term_map(goIds{i});
end
bestGoTerms = bestGoTerms';
bestGoIds = checkGos(sortInds)';
numIsesSorted = numUsedIsesAll(sortInds);
numGenesSorted = numDiffGenesAll(sortInds);


bestGoAllDataStruct = {};
for i = 1:length(goIds)
    bestGoTerms{i} = go_id_term_map(goIds{i});
    bestGoAllDataStruct{i, 1} = s(i);
    bestGoAllDataStruct{i, 2} = bestGoTerms{i};
    bestGoAllDataStruct{i, 3} = bestGoIds{i};
    
    bestGoAllDataStruct{i, 4} = numGenesSorted(i);
    bestGoAllDataStruct{i, 5} = numIsesSorted(i);
end

selectedGOs = [4984, -1, 0, 5667, 5216, 5244, 5249, 30165,61097,17146];%, 0, -1, , 5667, 55088, 6955, 45211, 5516];

tempEphysGeneCorrVec = [];
tempGONames = {};
for i = 1:length(selectedGOs)
    currGOTerm = selectedGOs(i);
    for j = 1:length(checkGos)
        if currGOTerm == checkGos{j}
            checkGOInd = j;
            break;
        end
    end
    tempEphysGeneCorrVec(i) = corrVals(j);
    if currGOTerm == 17146
        tempGONames{i} =  'NMDA glutamate receptor complex ';
    else
        tempGONames{i} = [' ' go_id_term_map(checkGos{j}) ' '];
    end
    tempCorrVecGeneN(i) = numDiffGenesAll(j);
end
figure;
bar(tempEphysGeneCorrVec, 'k');
set(gca, 'XTick', 1:length(tempEphysGeneCorrVec), 'XTickLabel', tempGONames);
box off;
hold on;
% for i = 1:length(tempEphysGeneCorrVec)
%     text(i-.1, tempEphysGeneCorrVec(i)+.03, num2str(tempCorrVecGeneN(i)));
% end
ylabel(' Electrophysiology - gene expression correlation (Pearson) ')
axis([0 length(selectedGOs)+1 0 max(tempEphysGeneCorrVec)*1.1])
xticklabel_rotate(1:length(tempEphysGeneCorrVec),90,tempGONames,'interpreter','none')
    
figure;
[c, x] = hist(corrVals(numDiffGenesAll>=3), 50, 'k');
bar(x, c, 'k');
ylabel('Count (Gene classes)')
box off;
xlabel(' Electrophysiology - gene expression correlation (Pearson) ')
% VGIon channels, neurotransmitter binding, bio_process, olfreceptors, nuc acid Trans factors, cytoskeleton, 

figure;
semilogx(numDiffGenesAll, corrVals, 'k.')
box off;
xlabel('Count (Genes per GO class)')
ylabel(' Electrophysiology - Gene expression correlation (Pearson) ')
hold on;
for i = 1:length(selectedGOs)
    text(tempCorrVecGeneN(i), tempEphysGeneCorrVec(i), tempGONames{i});
end

%% run cross validation to see which sets of linear model parameters are best
% e.g. aic vs bic and linear models or quadratic models
[allSses, fit_complexity_vec, fit_criteria] = metadata_fit_models_crossval(dataInputStruct);

metadata_cv_struct = [];
best_complexity = {};
for i = 1:6
    for j = 1:length(fit_complexity_vec)
        for k = 1:length(fit_criteria)
            metadata_cv_struct{i}.meanMat(j,k) = mean(squeeze(allSses(i,j,k,:)));
            metadata_cv_struct{i}.stdMat(j,k) = myStdErr(squeeze(allSses(i,j,k,:)));
        end
    end
    [~, tempInds] = min(metadata_cv_struct{i}.meanMat, [], 2);
    for k = 1:length(tempInds)
        tempVec(k) = metadata_cv_struct{i}.meanMat(k, tempInds(k));
    end
    [~,min_complexity_ind] = min(tempVec);
    [~,min_criteria_ind] = max(max(-metadata_cv_struct{i}.meanMat));
    best_complexity{i} = fit_complexity_vec{min_complexity_ind};
    best_criteria{i} = fit_criteria{min_criteria_ind};
end




% metadata_account_stepwise

%% plot results from cross-validated metadata

figure;
fit_complexity_vec = {' constant ',' linear ', ' purequadratic ', ' interactions ',' quadratic '};
fit_criteria = {'bic', 'aic'};
for i = 1:6
    subplot(2,3,i);

    errorbar(metadata_cv_struct{i}.meanMat(:,1), metadata_cv_struct{i}.stdMat(:,1));
    numEphysPts = sum(neuronUniqueArtMat(:,i));
    ylabel(' SSE (cross-validated)')
    title([ephysVarsNames{i} ' (n = ' num2str(numEphysPts) ') '])
    box off;
    set(gca, 'XTick', 1:5, 'XTickLabels', fit_complexity_vec)

    
end

%% run code for rerunning analysis on resampled data set - this takes a while too, make sure matlabpool on
ephys_analysis_resampled

% or just load data from ephys_neuron_clustering_resampling.mat

%% plot figures based on resampling analysis

figure;
for i = 1:numResponseVars
    subplot(2,3,i);
    boundedline(sampPctVec',meanEphysSampRsqNeuron(:,i), stdEphysSampRsqNeuron(:,i), 'k', 'alpha');
    hold on;
    boundedline(sampPctVec',meanEphysSampRsq(:,i), stdEphysSampRsq(:,i), 'r', 'alpha');   
    numEphysPts = sum(dataInputStructResamp(end).neuronUniqueArtMat(:,i));
    xlabel(' Resample fraction ');
    ylabel(' Model adj. R^2 ')
    title([ephysVarsNames{i} ' (n = ' num2str(numEphysPts) ') ']);
    axis([.18 1 0 1.02]);
end


barColors = ['k', 'r'];
figure; hold on;
barInd = 2;
barInds1 = [];
barInds2 = [];
for j = 1:6
    bar(barInd, neuronOnlyModelRsq(end,j), barColors(1));
    hold on;
%     errorbar(barInd, neuronOnlyModelRsq(end,j), stdEphysSampRsqNeuron(end-1,j), 'b');
    barInds1(end+1) = barInd;
    barInd = barInd + 1;
    bar(barInd, metadataModelRsq(end,j), barColors(2));
%     errorbar(barInd, metadataModelRsq(end,j), stdEphysSampRsq(end-1,j), 'b');
    barInds2(end+1) = barInd;
    barInd = barInd + 2;
    
end
errorbar(barInds1, neuronOnlyModelRsq(end,:), stdEphysSampRsqNeuron(end-1,:), '.b');
errorbar(barInds2, metadataModelRsq(end,:), stdEphysSampRsq(end-1,:), '.b');
ylabel(' Variance explained (adj. R^2) ')
axis([0 20 0 max(metadataModelRsq)*1.1]);
% set(gca, 'XTick', 2.5:3:18, 'XTickLabel', responseVarsNames)
format_ticks(gca,responseVarsNamesPlot, [], 2.5:3:18, [], [])
% set(gca, 'XTick', [])



% try to plot example distance matrices
simInds = [1 round(length(dataInputStructResamp)*.5) length(dataInputStructResamp)];
figure;
for i = 1:length(simInds)
    subplot(1,3,i);
    neuronInds = outStruct(simInds(i)).ephysNeuronInds;
    distMat = outStruct(simInds(i)).ephysDistMat;
    neuronInds = outStruct(simInds(i)).ephysNeuronInds;
    currNeuronNames = outStruct(simInds(i)).ephysNeuronNames;
    imagesc(distMat(neuronInds,neuronInds));
    colormap('redblue')
    caxis([0 max(max(ephysDistMat))]);
    set(gca, 'YTick', 1:length(keepData), 'YTickLabel', currNeuronNames(neuronInds));
    set(gca, 'XTick', []);
    title([' Resample fraction = ', num2str(dataInputStructResamp(simInds(i)).sampPct)]);
end


% 
figure; 
boundedline(sampPctVec, distMatCorrMeanVec, distMatCorrStdVec, 'k');   
hold on;
plot(sampPctVec, distMatCorrMeanVec, '.k');
xlabel(' Resample fraction ');
ylabel(' Electrophysiological distance similarity (Pearson corr) ');
axis([.18 1 0 1]);



%% integration with pvclust from R - need to run the analysis first in R
% how to integrate with R ward's distance:
recon = S.Recon;
keepData = recon;

n = neuronNames(neuronUseInds)';
r = recon';

distMat = squareform(pdist(keepData, metric));
figure;
imagesc(distMat(perm, perm));
colormap('redblue')

distVec = squareform(distMat);
% distVec = distVec.^(1/2);
linkageMethod = 'ward';
Z = linkage(distVec, linkageMethod);
% Z(:,3) = (Z(:,3)).^(1/2);
figure;
[H, t, perm] = dendrogram(Z,length(score), 'labels', neuronNames(neuronUseInds),   'Orientation', 'left');
axis([0 max(Z(:,3))+.1 -1 max(t)+1])

orig_neuron_tree =  phytree(Z, neuronNames(neuronUseInds));

% Pc = count ./ num_boots   % confidence probability (Pc)

[ptrs dist names] = get(orig_neuron_tree,'POINTERS','DISTANCES','NODENAMES');

num_seqs = length(neuronUseInds);
H = plot(orig_neuron_tree, 'BranchLabels', 'true')
for i = 1 : num_seqs -1  % for every branch
    branch_ptr = i + num_seqs;
%     names{branch_ptr} = [names{branch_ptr} ', confidence: ' num2str(100*Pc(i)) ' %'];
    %set(H.branchNodeLabels(i), 'String',num2str(round(Pc(i)*100)))
    set(H.branchNodeLabels(i), 'String',[num2str(round(au_pvalues(i)*100)),10,num2str(round(boot_prob(i)*100))]);
end

% neuronUseInds(phy_tree_inds)

phy_tree_names = get(orig_neuron_tree,'LEAFNAMES');
phy_tree_inds = [];
for i=1:length(phy_tree_names)
%      = strfind(neuronNames(neuronUseInds),phy_tree_names(i));
    phy_tree_inds(i)=find(ismember(neuronNames(neuronUseInds),phy_tree_names(i)));
end

figure; imagesc(keepData(phy_tree_inds, ephysUseInds));
colormap 'redblue'
colorbar;
set(gca, 'YTick', []);
set(gca, 'XTick', []);
% set(gca, 'YTick', 1:length(keepData), 'YTickLabel', phy_tree_names);
set(gca, 'XTick', 1:length(responseVarsNames), 'XTickLabel', responseVarsNames(ephysUseInds));
% box off;