function test_bug2225

% WALLTIME 00:10:00
% MEM 1gb

tic
for i=1:10000
  issue_warning
end % for
toc

end % main function

function issue_warning
ft_warning('this warning should not show too often');
end % subfunction
