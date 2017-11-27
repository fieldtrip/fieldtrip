function disp(obj)

% function disp(obj)
% Disp a configuration object. This function is generic, but it will be
% called also for derived objects except cfg_exbranch. It will first
% display fields inherited from cfg_item and then fields from the derived
% item.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: disp.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

sz = size(obj);
fprintf('%s object: ', class(obj));
if length(sz)>4,
    fprintf('%d-D\n',length(sz));
else
    for i=1:(length(sz)-1),
        fprintf('%d-by-',sz(i));
    end;
    fprintf('%d\n',sz(end));
end;
if prod(sz)==1,
    so = struct(obj);
    if isfield(so,'cfg_item')
        % derived objects have a field containing the base class
        disp(struct(so.cfg_item));
        % display additional fields, if any
        if numel(fieldnames(so)) > 1
            disp(rmfield(so,'cfg_item'));
        end;
    else
        % disp base object only
        disp(so);
    end;
end;
return;
