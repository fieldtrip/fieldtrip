function om_write_cond(condfile,c,names)
%   OM_WRITE_COND
%       [] = OM_WRITE_COND(CONDFILE,C)
%
%   Write conductivity file for OpenMEEG
%
%   Authors: Alexandre Gramfort alexandre.gramfort@inria.fr

ndomains = length(c);

if nargin < 3
    names = {};
    for k=1:ndomains
        names{k} = ['domain', num2str(k)];
    end
end

if ndomains ~= length(names)
    error('Number of conductivities is not equal to the number of domain names.');
end

cfid = fopen(condfile, 'w');
if cfid == -1
    error(['Failed to open file ''',condfile,'''.']);
end

fprintf(cfid,'# Properties Description 1.0 (Conductivities)\n\n');
fprintf(cfid,'Air         0.0\n');

for k=1:ndomains
    fprintf(cfid, '%s       %f\n', names{k}, c(ndomains-k+1));
end

fclose(cfid);

end %  function