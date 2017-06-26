function test_bug2225

% WALLTIME 00:10:00
% MEM 1gb

% This started failing on 22 June 2017 in relation to https://github.com/fieldtrip/fieldtrip/pull/445
% the "was_shown" output argument is not returned any more

counter=0;
for i=1:10
  was_shown=issue_warning();
  if was_shown
    counter=counter+1;
  end
end % for
if counter>1
  error('Too many repeated warnings shown');
end

end % main function

function was_shown=issue_warning
[unused,was_shown]=ft_warning('this warning should not show too often');
end % subfunction
