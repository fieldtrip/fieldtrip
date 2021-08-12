function test_ft_connectivity_powcorr_ortho

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_connectivity_powcorr_ortho

if nargout
  % assume that this is called by RUNTESTS
  tests = functiontests(localfunctions);
else
  % assume that this is called from the command line
  func = localfunctions;
  for i=1:numel(func)
    fprintf('evaluating %s\n', func2str(func{i}));
    feval(func{i});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testOptions(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nchan   = 10;
nrpttap = 100;

input = randn(nchan,nrpttap) + 1i*randn(nchan,nrpttap);

result = {};
result{end+1} = ft_connectivity_powcorr_ortho(input);
result{end+1} = ft_connectivity_powcorr_ortho(input, 'refindx', 1);
result{end+1} = ft_connectivity_powcorr_ortho(input, 'refindx', 1:2);
result{end+1} = ft_connectivity_powcorr_ortho(input, 'tapvec',  5*ones(1,20));
result{end+1} = ft_connectivity_powcorr_ortho(input, 'tapvec', 10*ones(1,10));
result{end+1} = ft_connectivity_powcorr_ortho(input, 'refindx', 1:2, 'tapvec',  5*ones(1,20));
result{end+1} = ft_connectivity_powcorr_ortho(input, 'refindx', 1:2, 'tapvec', 10*ones(1,10));

% all iterations were done with (slightly) different options, hence the results should not be equal
for i=1:numel(result)
  for j=(i+1):numel(result)
    assert(~isequaln(result{i}, result{j}), 'the results %d and %d should not be equal', i, j);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testImplementation(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mom1 = randn(1,100) + 1i*randn(1,100);
mom2 = randn(1,100) + 1i*randn(1,100);

c = ft_connectivity_powcorr_ortho([mom1; mom2], 'refindx', 1);
c = c(2,:);

% this is the computation according to the paper Nat Neuro 2012 Hipp et al.

% normalise the amplitudes
mom1n = mom1./abs(mom1);
mom2n = mom2./abs(mom2);
% rotate
mom12 = mom1.*conj(mom2n);
mom21 = mom2.*conj(mom1n);
% take the projection along the imaginary axis
mom12i = abs(imag(mom12));
mom21i = abs(imag(mom21));
% compute the correlation on the log transformed power values
c1 = corr(log10(mom12i.^2'), log10(abs(mom2).^2'));
c2 = corr(log10(mom21i.^2'), log10(abs(mom1).^2'));

% compare the FieldTrip implementation to that from the paper
assert(all(abs(c-(c1+c2)./2)<10*eps));
