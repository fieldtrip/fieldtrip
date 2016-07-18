function status = isequalfigure(a, b)

% ISEQUALFIGURE compares two figures and returns a boolean
%
% Use as
%  isequalfigure(a, b)
% with a and b being the figure handles.

fna = [tempname, '.png'];
fnb = [tempname, '.png'];

print(a, '-dpng', fna);
print(b, '-dpng', fnb);

mata = imread(fna);
matb = imread(fnb);

delete(fna);
delete(fnb);

status = isequal(mata, matb);

% the following is too sensitive
% md5a = ft_hash(fna);
% md5b = ft_hash(fnb);
% status = isequal(md5a, md5b);
