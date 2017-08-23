function test_xunit_stacktest

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_defaults ft_warning

ft_defaults

stack = dbstack('-completenames');
if size(stack) < 1
  fname = [];
  line = [];
  return;
end
i0 = 1;

fname = horzcat(stack(end).name);

fprintf('The fieldname would be tmpstrct.%s\n', fname);

tmpstrct = [];
setsubfield(tmpstrct, fname, [])
tmpstrct2 = [];
eval(['tmpstrct2.' fname ' =[]'])

for i=1:10
    warning(defaultId, 'Tthis is a test');
end
warning(defaultId, 'Tthis is another test');

warning(defaultId, '-clear');

fprintf('Completed!\n');

end
