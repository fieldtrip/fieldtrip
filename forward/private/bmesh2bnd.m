function out = bmesh2bnd(in)

% Helper to convert the fieldtrip bnd format for meshes into the bmesh
% format for the HBF toolbox, and vice versa.
%
% George O'Neill

n = numel(in);

switch class(in)
    case 'cell'
        % convert bmeshes to bnds
       out = cellfun(@(c) struct('tri', c.e, 'pos', c.p),...
           in, 'UniformOutput', true);
    case 'struct'
        % convert bmeshes to bmeshes
        out = arrayfun(@(s) struct('e', s.tri, 'p', s.pos),...
            in, 'UniformOutput', false);
end

