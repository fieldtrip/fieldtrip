function [a, b] = isintent(this,intent)
% Correspondance between fieldnames and NIfTI intents
% FORMAT ind = isintent(this,intent)
% this    -  GIfTI object
% intent  -  fieldnames
% a       -  indices of found intent(s)
% b       -  indices of dataarrays of found intent(s)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id$

a = [];
b = [];
if ischar(intent), intent = cellstr(intent); end
c = cdata;
for i=1:length(this(1).data)
    switch this(1).data{i}.attributes.Intent(14:end)
        case 'POINTSET'
            [tf, loc] = ismember('vertices',intent);
            if tf
                a(end+1) = loc;
                b(end+1) = i;
            end
            [tf, loc] = ismember('mat',intent);
            if tf
                a(end+1) = loc;
                b(end+1) = i;
            end
        case 'TRIANGLE'
            [tf, loc] = ismember('faces',intent);
            if tf
                a(end+1) = loc;
                b(end+1) = i;
            end
        case 'VECTOR'
            [tf, loc] = ismember('normals',intent);
            if tf
                a(end+1) = loc;
                b(end+1) = i;
            end
        case {'NONE', 'LABEL', 'SHAPE', 'TIME_SERIES', 'RGB_VECTOR', ...
                'RGBA_VECTOR' c{:}}
            [tf, loc] = ismember('cdata',intent);
            if tf
                a(end+1) = loc;
                b(end+1) = i;
            end
        otherwise
            fprintf('Intent %s is ignored.\n',this.data{i}.attributes.Intent);
    end
end
[d,i] = unique(a,'first');
if length(d) < length(a)
    %warning('Several fields match intent type. Using first.');
    a = a(i);
    %b = b(i);
end

function c = cdata

c = {
'CORREL'
'TTEST'
'FTEST'
'ZSCORE'
'CHISQ'
'BETA'
'BINOM'
'GAMMA'
'POISSON'
'NORMAL'
'FTEST_NONC'
'CHISQ_NONC'
'LOGISTIC'
'LAPLACE'
'UNIFORM'
'TTEST_NONC'
'WEIBULL'
'CHI'
'INVGAUSS'
'EXTVAL'
'PVAL'
'LOGPVAL'
'LOG10PVAL'
'ESTIMATE'
'LABEL'
'NEURONAMES'
};
