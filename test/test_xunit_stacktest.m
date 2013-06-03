function test_xunit_stacktest


% TEST test_xunit_stacktest
% TEST ft_defaults warning_once

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

fprintf('Completed!\n');

end