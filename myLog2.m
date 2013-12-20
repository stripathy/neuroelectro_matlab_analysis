function [output] = myLog2(input)

output = log2(input);
inds = isinf(output);
output(inds) = 0;

inds = isnan(output);
output(inds) = 0;