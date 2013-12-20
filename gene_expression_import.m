
fName = 'data/brain_region_info.csv';
fid = fopen(fName,'r');  % Open text file

[ndata, text_data, alldata] = xlsread(fName);

% predictorVarsNames = {'NeuronType';	'Species';	'Strain';	'ElectrodeType'; 'PrepType'; ...
%     'JxnPotential'; 'JxnOffset'; 'Age'; 'Temp'; 'Weight'; 'PubYear'};
% responseVarsNames = {'rmp'; 'ir'; 'tau'; 'amp'; 'hw'; 'thresh'};
% 
% contVars = {'rmp'; 'ir'; 'tau'; 'amp'; 'hw'; 'thresh'; 'JxnOffset'; 'Age'; 'Temp'; 'Weight'; 'PubYear'};
% 
% catVarInds = [1 1 1 1 1 0 0 0];
headers = alldata(1,:);
regionStruct = {};
for i = 2:size(alldata,1)
    for j = [1 2 3 4 5]
        if j == 1
            regionStruct{i-1}.regionid = alldata{i,j};
        elseif j == 2
            regionStruct{i-1}.acronym = alldata{i,j};
        elseif j == 3
            regionStruct{i-1}.name = alldata{i,j};
        elseif j == 4
            regionStruct{i-1}.treedepth = alldata{i,j};
        elseif j == 5
            regionStruct{i-1}.color = alldata{i,j};
        end
    end
end

fName = 'data/image_series_info.csv';
fid = fopen(fName,'r');  % Open text file

[ndata, text_data, alldata] = xlsread(fName);

headers = alldata(1,:);
iseStruct = {};
for i = 2:size(alldata,1)
    for j = [1 2 3 4 5]
        if j == 1
            iseStruct{i-1}.allen_id = alldata{i,j};
        elseif j == 2
            iseStruct{i-1}.acronym = alldata{i,j};
        elseif j == 3
            iseStruct{i-1}.gene_name = alldata{i,j};
        elseif j == 4
            iseStruct{i-1}.ise_id = alldata{i,j};
        elseif j == 5
            iseStruct{i-1}.plane = alldata{i,j};
        elseif j == 6
            iseStruct{i-1}.entrez_id = alldata{i,j};
        elseif j == 7
            data = alldata{i,j};
            if strcmp('FALSE', data)
                data = 0;
            else
                data = 1;
            end
            iseStruct{i-1}.valid = data;
        end
    end
end

planeVec = zeros(length(iseStruct), 1);
for i = 1:length(planeVec)
    if strcmp(iseStruct{i}.plane, 'sagittal')
        planeVec(i) = 2;
    else
        planeVec(i) = 1;
    end
end


regionNames= {};
for i = 1:length(regionStruct)
    regionNames{i} = regionStruct{i}.name;
end


fName = 'data/neuron_region_assignment_modified.csv';
fid = fopen(fName,'r');  % Open text file

[ndata, text_data, alldata] = xlsread(fName);

clear text 
headers = alldata(1,:);
neuronStruct = {};
for i = 2:size(alldata,1)
    neuronStruct{i-1}.name = [alldata{i,1} ' '];
    regionListNames = {};
    regionListIds = [];
    cnt = 1;
    for j = 4:size(alldata,2)
        if ~isnan(alldata{i,j})
            regionListNames{cnt} = alldata{i,j};
            regionId = find(strcmp(alldata{i,j}, regionNames));
            regionListIds(cnt) = regionId;
            cnt = cnt + 1;
        end
    end
    if cnt > 1
        restrictInd = alldata{i,2};
        permissiveInd =  alldata{i,3};
        neuronStruct{i-1}.regionListNames = regionListNames;
        if restrictInd == 1
            neuronStruct{i-1}.regionListIdsRest = regionListIds;
        end
        if permissiveInd == 1
            neuronStruct{i-1}.regionListIdsPerm = regionListIds;
        end
    end
end

fName = 'data/gene_energy_mat.csv';
fid = fopen(fName,'r');  % Open text file
[geneEnergyMat] = csvread(fName);


gene_ise_index_map = containers.Map('KeyType', 'char', 'ValueType','any');
for i = 1:length(iseStruct)
    currIse = iseStruct{i};
    gene_name = char(currIse.acronym);
    if isKey(gene_ise_index_map,gene_name)
        gene_ise_index_map(gene_name) = [gene_ise_index_map(gene_name) i];
    else
        gene_ise_index_map(gene_name) = i;
    end
end

