
dataUseInds = 1:length(NeuronType);
[~, dataUseInds] = unique(Title, 'stable');

responseVars = [rmp	ir tau	amp	hw thresh];
predictorVarsCat = [NeuronType	Species	Strain	ElectrodeType PrepType JxnPotential];
% predictorVarsCont = [Age Temp Weight];
predictorVarsCont = [Age Temp];

predictorVars = []
for i = 1:size(predictorVarsCat, 2)
    predictorVars(:,end+1) = toCategorical(predictorVarsCat(dataUseInds,i));
end

for i = 1:size(predictorVarsCont, 2)
    predictorVars(:,end+1) = predictorVarsCont(dataUseInds, i);
end
% predictorVarsNames = {'NeuronType';	'Species';	'Strain';	'ElectrodeType'; 'PrepType'; 'Age'; 'Temp'; 'Weight'};
predictorVarsNames = {'NeuronType';	'Species';	'Strain';	'ElectrodeType'; 'JxnPotential'; 'PrepType'; 'Age'; 'Temp'};
responseVarsNames = {'rmp'; 'ir'; 'tau'; 'amp'; 'hw'; 'thresh'};
responseVarsNamesPlot = {' V_{rest} '; ' R_{input} '; ' \tau_{m} '; ' AP_{amp} '; ' AP_{hw} '; ' AP_{thr} '};

% catVarInds = [1 1 1 1 1 0 0 0];
catVarInds = [1 1 1 1 1 1 0 0];

figure;
% species/strain
metadataVals = Strain;
metadataInd = 3;
numMetadataVals = length(unique(metadataVals));
metadataXVals = unique(metadataVals, 'stable');

orderInds = [3 2 5 7 1 4 6];

subplot(2,5, 1:2);

[count, binCenters] = hist(predictorVars(:, metadataInd), numMetadataVals);
bar(binCenters, count(orderInds), 'k');
set(gca, 'XTick', binCenters, 'XTickLabel', metadataXVals(orderInds));
ylabel('Count (articles) ');
box off;
axis([1 nnz(count>0) 0 max(count)*1.05])

% Electrode type
% orderInds = [1 2 3];
metadataVals = ElectrodeType;
metadataInd = 4;
numMetadataVals = length(unique(metadataVals));
metadataXVals = unique(metadataVals, 'stable');

subplot(2,5, 3);

[count, binCenters] = hist(predictorVars(:, metadataInd), numMetadataVals);
bar(1:nnz(count>0), count(count>0), 'k');
set(gca, 'XTick', 1:nnz(count>0), 'XTickLabel', metadataXVals(1:nnz(count>0)));
% set(gca, 'XTick', binCenters, 'XTickLabel', metadataXVals);
box off;
axis([1 nnz(count>0) 0 280])

% PrepType
metadataVals = PrepType;
metadataInd = 5;
numMetadataVals = length(unique(metadataVals));
metadataXVals = unique(metadataVals, 'stable');
% metadataXVals = {'In vitro '; 'In vivo '; 'Culture ';};

subplot(2,5, 4);

[count, binCenters] = hist(predictorVars(:, metadataInd), numMetadataVals);
bar(1:nnz(count>0), count(count>0), 'k');
set(gca, 'XTick', 1:nnz(count>0), 'XTickLabel', metadataXVals(1:nnz(count>0)));
% set(gca, 'XTick', binCenters, 'XTickLabel', metadataXVals);
box off;
% axis([1 nnz(count>0) 0 280])
axis tight


% Junction Potential
orderInds = [1 2 3];
metadataVals = JxnPotential;
metadataInd = 6;
numMetadataVals = length(unique(metadataVals));
metadataXVals = unique(metadataVals, 'stable');
subplot(2,5, 5);

[count, binCenters] = hist(predictorVars(:, metadataInd), numMetadataVals);
bar(binCenters, count(orderInds), 'k');
set(gca, 'XTick', binCenters, 'XTickLabel', metadataXVals(orderInds));
box off;
axis([1 nnz(count>0) 0 280])

% Age
metadataVals = Age(Age<300);
metadataInd = 7;
numMetadataVals = 25;

subplot(2,5, 6:7);

[count, binCenters] = hist(metadataVals, numMetadataVals);
bar(binCenters, count, 'k');
% set(gca, 'XTick', binCenters, 'XTickLabel', metadataXVals);
xlabel('Animal age (days) ');
ylabel('Count (articles) ');
box off;
axis([min(binCenters)-1 250 0 max(count)*1.1])

% Temp
metadataVals = Temp;
metadataInd = 8;
numMetadataVals = 17;

subplot(2,5, 9:10);

[count, binCenters] = hist(predictorVars(:, metadataInd), numMetadataVals);
bar(binCenters, count, 'k');
% set(gca, 'XTick', binCenters, 'XTickLabel', metadataXVals);
ylabel('Count (articles) ');
xlabel('Recording temperature (°C) ');
box off;
axis([min(binCenters)-1 max(binCenters)+1 0 max(count)*1.1])

% % Neuron Type
% metadataVals = NeuronType;
% metadataInd = 1;
% numMetadataVals = 10;
% neuronFreqs = [];
% for i = 1:length(unique(predictorVars(:, metadataInd)))
%     neuronFreqs(i) = sum(predictorVars(:, metadataInd) == i);
% end

% figure;
% [count, binCenters] = hist(neuronFreqs, max(unique(neuronFreqs)));
% bar(binCenters, count, 'k');
% ylabel('Count (articles) ');
% xlabel('Count (neuron types) ');
% box off;
%%
% Neuron Type
metadataVals = NeuronType;
metadataInd = 1;
numMetadataVals = 10;
neuronFreqs = [];
for i = 1:length(unique(predictorVars(:, metadataInd)))
    neuronFreqs(i) = sum(predictorVars(:, metadataInd) == i);
end

figure;
subplot(1,2,1);
[count, binCenters] = hist(neuronFreqs, max(unique(neuronFreqs)));
bar(binCenters, count, 'k');
xlabel('Count (articles) ');
ylabel('Count (neuron types) ');
box off;

% EphysHisto

subplot(1,2,2);
ephysArticleCounts = sum(~isnan(responseVars));
bar(1:size(responseVars,2), ephysArticleCounts, 'k');
ylabel('Count (measurements) ');
format_ticks(gca,responseVarsNamesPlot, [], [], [], 90)
set(gca, 'XTick', [])

% set(gca, 'XTick', 1:size(responseVars,2), 'XTickLabel', responseVarsNames)

axis([0 size(responseVars,2)+1 0 max(ephysArticleCounts)*1.1])
box off;

%%
figure;
predictorVarsNames = {'NeuronType';	'Species';	'Strain';	'ElectrodeType'; 'JxnPotential'; 'PrepType'; 'Age'; 'Temp'};
acc = [168/276, 141/249, 148/273, 74/94, 156/278, 27/246, 89/251];
bar(acc*100, 'k');
format_ticks(gca,predictorVarsNames(2:end), [], [], [], 90)
ylabel('Assignment accuracy (%) ');
box off;


