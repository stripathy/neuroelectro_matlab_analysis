% loads in neuroelectro data from spreadsheets and puts into variables and
% format

fName = 'data/article_ephys_metadata_summary_10_24_13.csv';
fid = fopen(fName,'r');  % Open text file

[ndata, text, alldata] = xlsread(fName);

predictorVarsNames = {'NeuronType';	'Species';	'Strain';	'ElectrodeType'; 'PrepType'; ...
    'JxnPotential'; 'JxnOffset'; 'Age'; 'Temp'; 'Weight'; 'PubYear'};
responseVarsNames = {'rmp'; 'ir'; 'tau'; 'amp'; 'hw'; 'thresh'};

contVars = {'rmp'; 'ir'; 'tau'; 'amp'; 'hw'; 'thresh'; 'JxnOffset'; 'Age'; 'Temp'; 'Weight'; 'PubYear'};

% variables indicating categorical variables
catVarInds = [1 1 1 1 1 0 0 0];

headers = alldata(1,:);
for i = 1:length(headers)
    currHeader = headers{i};
    currVec = alldata(:,i);
    if sum(strcmp(contVars, currHeader)) > 0
        
        dataVec = [];
        for j = 2:length(currVec)
            if strcmp('NaN', currVec{j})
                dataVec(j-1,1) = NaN;
            elseif strcmp('nan', currVec{j})
                dataVec(j-1,1) = NaN;
            else
                dataVec(j-1,1) = double(currVec{j});
            end
        end
        eval([currHeader ' = dataVec;']);
    else
        strVec = {};
        for j = 2:length(currVec)
            if isnan(currVec{j})
                strVec{j-1,1} = '';
            else
                strVec{j-1,1} = [currVec{j} ' '];
            end
        end
        eval([currHeader ' = strVec;']);
    end
end
% ir = log10(ir);
% tau = log10(tau);
% Age = log10(Age);
for i = 1:length(JxnPotential)
    if strcmp('', JxnPotential(i))
        JxnPotential{i} = 'Not reported ';
    end
end

clear text

for i = 1:length(NeuronType)
    if strcmp('Nucleus accumbens shell neuron ', NeuronType(i))
        NeuronType{i} = 'Nucleus accumbens medium spiny neuron ';
    elseif strcmp('Nucleus accumbens core neuron ', NeuronType(i))
        NeuronType{i} = 'Nucleus accumbens medium spiny neuron ';
    elseif strcmp('Neocortex pyramidal cell ', NeuronType(i))
        NeuronType{i} = 'Neocortex uncharacterized cell ';
    end
end

% gives guinea pigs the strain guinea pigs
for i = 1:length(Species)
    if strcmp('Guinea Pigs ', Species{i})
        Strain{i} = 'Guinea Pigs ';
    end
end

for i = 1:length(Strain)
    if strcmp('', Strain{i})
        Strain{i} = 'Other ';
    end
end
        
    
