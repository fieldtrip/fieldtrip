function [a, b] = isintent(this,intent)
% Correspondance between fieldnames and NIfTI intent codes
% FORMAT ind = isintent(this,intent)
% this    -  GIfTI object
% intent  -  fieldnames
% a       -  indices of found intent(s)
% b       -  indices of dataarrays of found intent(s)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: isintent.m 6345 2015-02-20 12:25:50Z guillaume $

a = [];
b = [];
if ischar(intent), intent = cellstr(intent); end
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
        case 'NODE_INDEX'
            [tf, loc] = ismember('indices',intent);
            if tf
                a(end+1) = loc;
                b(end+1) = i;
            end
        case cdata
            [tf, loc] = ismember('cdata',intent);
            if tf
                a(end+1) = loc;
                b(end+1) = i;
            end
            if strcmp(this(1).data{i}.attributes.Intent(14:end),'LABEL')
                [tf, loc] = ismember('labels',intent);
                if tf
                    a(end+1) = loc;
                    b(end+1) = i;
                end
            end
        otherwise
            fprintf('Intent %s is ignored.\n',this.data{i}.attributes.Intent);
    end
end
%[d,i] = unique(a);
%if length(d) < length(a)
%    warning('Several fields match intent type. Using first.');
%    a = a(i);
%    b = b(i);
%end

function c = cdata

c = {
'NONE'
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
'GENMATRIX'
'SYMMATRIX'
'DISPVECT'
'QUATERNION'
'DIMLESS'
'TIME_SERIES'
'RGB_VECTOR'
'RGBA_VECTOR'
'SHAPE'
'CONNECTIVITY_DENSE'
'CONNECTIVITY_DENSE_TIME'
'CONNECTIVITY_PARCELLATED'
'CONNECTIVITY_PARCELLATED_TIME'
'CONNECTIVITY_CONNECTIVITY_TRAJECTORY'
};
