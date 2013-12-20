


numValidNeurons = length(validNeuronInds);


if go_id < 1
    go_name = 'all genes';
else
    go_name = GO(go_id).term.name;
end

coronal_only = 0;
average_ise_genes = 1;
minExpThreshVoxel = 2^-Inf;
minExpThreshExpt = 2^-Inf;

useGeneNames = go_gene_map_full(go_id);
geneSimInds = [];
geneUseInds = [];
for i = 1:length(useGeneNames)
    if iscellstr(useGeneNames(i))
        gene_name = char(useGeneNames(i));
    
        if isKey(gene_ise_index_map, gene_name)
            currIseInds = gene_ise_index_map(gene_name);
            for j = 1:length(currIseInds)
                meanExptEnergy = geneEnergyMat(currIseInds(j), 1);
                if meanExptEnergy >= minExpThreshExpt
                    geneUseInds(end+1) = currIseInds(j);
                    geneSimInds(end+1) = i;
                end
            end
        end
    end
end

if coronal_only == 1
    geneUseInds = intersect(geneUseInds, find(planeVec == 1));
end
% geneUseInds = 1:26065;


currGeneMat = zeros(numValidNeurons, length(geneUseInds));

for i = 1:numValidNeurons
    regionInds = neuronStruct{validNeuronInds(i)}.regionListIdsPerm;
    tempGeneMat = geneEnergyMat(geneUseInds, regionInds);
    tempGeneMat(tempGeneMat <= minExpThreshVoxel) = minExpThreshVoxel;
    tempGeneMat = nanmean(myLog2(tempGeneMat),2);
    currGeneMat(i,:) = tempGeneMat';
end

uniqueGeneInds = unique(geneSimInds, 'stable');
numUniqueGenes = length(uniqueGeneInds);

if average_ise_genes == 1
    avgGeneMat = zeros(numValidNeurons, numUniqueGenes);
    for i = 1:numUniqueGenes
        avgGeneMat(:,i) = nanmean(currGeneMat(:,geneSimInds==uniqueGeneInds(i)), 2);
    end
else
    avgGeneMat = currGeneMat;
end

currGeneMat = avgGeneMat;
% currGeneMat = zscore(avgGeneMat);




keepData = currGeneMat;

metric = 'euclidean';
if strcmp(metric, 'euclidean')
    linkageMethod = 'average';
else
    linkageMethod = 'average';
end



% cg = clustergram(keepData, 'RowLabels', neuronNames(validNeuronInds), ...
%     'Cluster', 3, 'Standardize', 'none', 'OptimalLeafOrder', 'false', 'Colormap', 'redblue');

% Z = linkage(keepData, linkageMethod, metric);
% figure;
% [H, t, perm] = dendrogram(Z,numValidNeurons);
% xticklabel_rotate

% figure; imagesc(keepData(perm, :));
% set(gca, 'YTick', 1:length(keepData), 'YTickLabel', neuronNames(m(perm)));
% set(gca, 'XTick', 1:length(responseVarsNames(ephysUseInds)), 'XTickLabel', responseVarsNames(ephysUseInds));
% % xticklabel_rotate
% colorbar;

geneDistMat = squareform(pdist(keepData, metric));
% figure;
% imagesc(distMat(perm, perm));
% colorbar;
% set(gca, 'YTick', 1:length(keepData), 'YTickLabel', neuronNames(validNeuronInds(perm)));
% set(gca, 'XTick', 1:length(keepData), 'XTickLabel', neuronNames(m(perm)));

geneNeuronNames = neuronNames(validNeuronInds);
geneRegionNames = regionNames(validNeuronInds);


%%

allNeuronNames = ephysNeuronNames;
sharedNeuronTypes = {};
ephysSharedInds = [];
genesSharedInds = [];
for i = 1:length(allNeuronNames)
    if  sum(strcmp(allNeuronNames(i), geneNeuronNames')) > 0 && sum(strcmp(allNeuronNames(i), ephysNeuronNames)) > 0
        genesSharedInds(end+1) = find(strcmp(allNeuronNames(i), geneNeuronNames'));
        ephysSharedInds(end+1) = find(strcmp(allNeuronNames(i), ephysNeuronNames'));
        sharedNeuronTypes(end+1) = allNeuronNames(i);
    end
end

ephysSubDistMat = ephysDistMat(ephysSharedInds,ephysSharedInds);
geneSubDistMat = geneDistMat(genesSharedInds,genesSharedInds);

% Z = linkage(ephysSubDistMat, linkageMethod, metric);
% [H, t, perm] = dendrogram(Z,length(sharedNeuronTypes));

figure;
imagesc(ephysSubDistMat);
set(gca, 'YTick', 1:length(sharedNeuronTypes), 'YTickLabel', sharedNeuronTypes);
set(gca, 'XTick', []);
h = colorbar('SouthOutside');
xlabel(h, ' Electrophysiological distance (Euclidean) ')
colormap redblue;

figure;
imagesc(geneSubDistMat);
set(gca, 'YTick', 1:length(sharedNeuronTypes), 'YTickLabel', geneRegionNames(genesSharedInds), 'YAxisLocation', 'right');
set(gca, 'XTick', []);
h = colorbar('SouthOutside');
xlabel(h, ' Gene expresion distance (Euclidean) ')
colormap redblue;
title(['GO Term: ', go_name]);
% caxis([.6 1.2])

ephysPairwiseDistVals = squareform(ephysSubDistMat)';
genePairwiseDistVals = squareform(geneSubDistMat)';



figure;
plot(ephysPairwiseDistVals, genePairwiseDistVals, 'k.')
xlabel(' Electrophysiological distance ')
ylabel(' Gene expression distance ')


% fitInds = (ephysPairwiseDistVals < 2);
% figure; plot(ephysPairwiseDistVals(fitInds), genePairwiseDistVals(fitInds), '.');
% corr(ephysPairwiseDistVals(fitInds), genePairwiseDistVals(fitInds))

% [r, p] = corr(squareform(ephysSubDistMat)', squareform(geneSubDistMat)', 'type', 'Spearman')

[r,p] = f_mantel(ephysSubDistMat, geneSubDistMat,1,mantel_resamples)

    



        