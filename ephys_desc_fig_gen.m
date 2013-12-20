% load ephys + metadata from spreadsheet
% ephys_metadata_import_fxn
ephysData = [rmp	ir tau	amp	hw thresh];
% plot simple summary data for neurons
useNeuronNames = {'Cerebellum Purkinje cell ', 'Hippocampus CA1 pyramidal cell ', 'Neocortex basket cell ', ...
    'Ventral tegmental area dopamine neuron ', 'Neostriatum medium spiny neuron '};
useNeuronShortNames = {' Purk ', ' CA1, pyr ', ' CTX, basket ', ' VTA DA ', ' STR, MSN '};
neuronColors = {'r', 'b', 'm', 'c', 'g'};
numUseNeurons = length(useNeuronNames);
ephysVarsNamesShort = {'rmp'; 'ir'; 'tau'; 'amp'; 'hw'; 'thresh'};
ephysVarsNames = {' Resting membrane potential'; ' Input resistance'; ' Membrane time constant';...
    ' Spike amplitude'; ' Spike half-width'; ' Spike threshold'};
ephysVarsUnits = {'mV'; 'MOhms'; 'ms';...
    'mV'; 'ms'; 'mV'};

% data point filtering
dataUseInds = log10(Age) > .67;
[neuronCategorical, neuronNames] = grp2idx(NeuronType(dataUseInds));

ephysInds = 1:6;
figure;
for i = 1:length(ephysInds)
    subplot(2, 3, i); hold on;
    for j = 1:numUseNeurons
        tempNeuronInds = strcmp(useNeuronNames(j), NeuronType);
%         tempNeuronInds = tempNeuronInds(dataUseInds);
        jitterInds1 = repmat(j, sum(tempNeuronInds), 1) + randn(sum(tempNeuronInds), 1)*.1;
        if i ~= 2
            plot(jitterInds1, ephysData(tempNeuronInds, i), [neuronColors{j} 'o']);
        else
            semilogy(jitterInds1, ephysData(tempNeuronInds, i), [neuronColors{j} 'o']);
        end
    end
    set(gca, 'XTick', [1:numUseNeurons], 'XTickLabel', useNeuronShortNames);
%     xticklabel_rotate;
    ylabel([ephysVarsNames{i} ' (' ephysVarsUnits{i} ') ']);
end