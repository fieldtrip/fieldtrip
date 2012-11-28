function test_suite = test_bug1618
% TEST test_bug1618
%

% add xunit to path
ft_hastoolbox('xunit',1);
initTestSuite;  % for xUnit

function test_no_nan
  % TODO: load example data, test for absence of NAN.
  h = ft_read_header('data_bug1618/bug1618.dat');
  h
  h.chanunit
  h.chantype

  X = ft_read_data('data_bug1618/bug1618.dat');
  assert(~any(isnan(X(:))), 'NaN values are not expected for this dataset!');

