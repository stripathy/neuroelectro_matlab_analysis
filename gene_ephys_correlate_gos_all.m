minIseExpts = 10;
minUniqueGenes = 10;

% minExpThreshVoxel = 2^-Inf;
% minExpThreshExpt = 2^-Inf;
% average_ise_genes = 1;

average_ise_genes = 1;
minExpThreshVoxel = 2^-Inf;
minExpThreshExpt = 2^-Inf;

% go_id = 51539;
checkGos = keys(go_gene_map_full);
corrVals = zeros(length(checkGos), 1);
corrValsSpearman = zeros(length(checkGos), 1);
corrValsPearson = zeros(length(checkGos), 1);
numDiffGenesAll = zeros(length(checkGos), 1);
numUsedIsesAll = zeros(length(checkGos), 1);
parfor t = 1:length(checkGos)
    
    go_id = checkGos{t};
    coronal_only = 0;
    
    useGeneNames = go_gene_map_full(go_id);

    geneUseInds = [];
    geneSimInds = [];
    for i = 1:length(useGeneNames)
        if iscellstr(useGeneNames(i)) % check that the gene name is valid
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
    
    if length(geneUseInds) < minIseExpts
        continue
    end
    
    if coronal_only == 1
        geneUseInds = intersect(geneUseInds, find(planeVec == 1));
    end
    % geneUseInds = 1:26065;
    
    numUsedIses = length(geneUseInds);
    
    
    
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
    
    if numUniqueGenes < minUniqueGenes
        continue
    end
    
    
    avgGeneMat = zeros(numValidNeurons, numUniqueGenes);
    if average_ise_genes == 1 % average over experiments corresponding to same gene
        if numUniqueGenes < 3
            continue;
        end
        for i = 1:numUniqueGenes
            avgGeneMat(:,i) = nanmean(currGeneMat(:,geneSimInds==uniqueGeneInds(i)), 2);
        end
    end

%     currGeneMat = zscore(avgGeneMat);% + randn(size(avgGeneMat))*.001;
    currGeneMat =avgGeneMat;
    
    keepData = currGeneMat;
    
    metric = 'euclidean';
    if strcmp(metric, 'euclidean')
        linkageMethod = 'average';
    else
        linkageMethod = 'average';
    end
    
    distMat = squareform(pdist(keepData, metric));
    geneDistMat = distMat;
    geneNeuronNames = neuronNames(validNeuronInds);
    
    %%
    allNeuronNames = union(geneNeuronNames, ephysNeuronNames);
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
    
    
%     [r,p] = f_mantel(ephysSubDistMat(perm, perm), geneSubDistMat(perm, perm),1,100);
    [r,p] = f_mantel(ephysSubDistMat, geneSubDistMat,1,100);
    [r_spearman] = corr(squareform(ephysSubDistMat)', squareform(geneSubDistMat)', 'type', 'Spearman');
    [r_correlation] = corr(squareform(ephysSubDistMat)', squareform(geneSubDistMat)', 'type', 'Pearson');
    corrVals(t) = r;
    corrValsSpearman(t) = r_spearman;
    corrValsPearson(t) = r_correlation;
    numUsedIsesAll(t) = numUsedIses
    numDiffGenesAll(t) = numUniqueGenes;
    
end


