function [label_num, label_names, Nclass] = label2num(label, label_names)
% Create label vector consisting of the number from labels.
%  
% Note that elements of output labels starts from 1. 
%
% -- Usage
% [label_num, label_names, Nclass] = label2num(label, label_names)
%
% -- Example 
% > label = {'red', 'green', 'green', 'red'}
% > label_num = label2num(label)
% > [1; 2; 2; 1];
%
% -- Input
% label : label vector 
% 
% 2009/08/10 OY 
% * add the second input argument as an optional input
% If "label_names" is provided, "label" is transrated to "label_num" 
% according to this "label_names". The order of "label_names" 
% is important. 
% ex. label_names = {'red','green'}
%     label = {'green','green','green'} --> label_num = [2 2 2]
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.


if (ischar(label))
   label = cellstr(label);
end
if (size(label, 1) == 1)
   label = label';
end

%
if nargin < 2
label_names = unique(label);
end
Nclass = length(label_names);
Nsamp = length(label);

label_num = NaN * ones(Nsamp,1);

for ii = 1 : Nclass
    if iscell(label_names)
        ix = find(label == label_names{ii});
    else
        ix = find(label == label_names(ii));
    end
    
    label_num(ix) = ii;
    
end
    

