function [output] = toCategorical(input)

[output]=grp2idx(input);

for i = 1:length(input)
    if cellfun('isempty', input(i))
        output(i) = NaN;
    end
end
