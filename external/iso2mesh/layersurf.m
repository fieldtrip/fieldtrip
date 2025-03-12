function [face, labels] = layersurf(elem, varargin)
%
% face=layersurf(elem, opt)
%    or
% [face,labels]=layersurf(elem,'option1',value1,'option2',value2,...)
%
% process a multi-layered tetrahedral mesh, like a brain mesh, to extract
% the layer surface meshes with one enclosed by another
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% input:
%    elem: an Nx5 integer array, representing the tetrahedral mesh element
%         list. The first 4 columns represent the tetrahedral element node
%         indices; the last column represents tissue labels.
%    opt: (optional) a struct or pairs of names/values to provide
%         additional options; accepted options include
%         'order': ['>='] or '=', '<=' . if set to '>=' (default), the
%             outmost layer has the lowest label count; if '<=', innermost
%             is lowest; if '=', surface of each label is extracted
%             individually - meaning that inner surfaces will have two
%             duplicates
%         'innermost': [0] or an array of labels. The labels defined in
%             this list will be treated as the innermost regions, and its
%             boundary will be extracted using the '==' test (order='=').
%             by default; label 0 is assumed to be innermost (i.e. nothing
%             is enclosed inside).
%         'unique': [0] or 1. if 1, remove duplicated triangles; if 0, keep
%             all triangles.
%         'occurence': ['first'] or 'last'. if 'first', the unique operator
%             keeps a duplicated triangle with the lowest label number;
%             otherwise, a repeated triange keeps the highest label number.
%
% output:
%    vol: a 2-D or 3-D array of the same type/size of the input arrays. The
%         label for each voxel is determined by the index to the highest
%         value in TPM of the same voxel. If a voxel is a background voxel
%         - i.e. zeros for all TPMs, it stays 0
%    names: a cell array storing the names of the labels (if input is a
%         struct), the first string is the name for label 1, and so on
%
% -- this function is part of brain2mesh toolbox (http://mcx.space/brain2mesh)
%    License: GPL v3 or later, see LICENSE.txt for details
%

opt = varargin2struct(varargin{:});
outsideislower = jsonopt('order', '>=', opt);
dounique = jsonopt('unique', false, opt);
innermost = jsonopt('innermost', 0, opt);
occurence = jsonopt('occurence', 'first', opt);

labels = sort(unique(elem(:, 5)));
face = [];
for i = 1:length(labels)
    if (strcmp(outsideislower, '>=') && ~ismember(labels(i), innermost))
        newface = volface(elem(elem(:, 5) >= labels(i), 1:4));
    elseif (strcmp(outsideislower, '<=') && ~ismember(labels(i), innermost))
        newface = volface(elem(elem(:, 5) <= labels(i), 1:4));
    else
        newface = volface(elem(elem(:, 5) == labels(i), 1:4));
    end
    newface(:, 4) = labels(i);
    face = [face; newface];
end

if (dounique)
    face(:, 1:3) = sort(face(:, 1:3), 2);
    [uniqface, idx] = unique(face(:, 1:3), 'rows', occurence);
    face = [uniqface, face(idx, end)];
end
