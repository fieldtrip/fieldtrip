function [indexed, probabilistic] = determine_segmentationstyle(segmentation, fn, dim)

% DETERMINE_SEGMENTATIONSTYLE is a helper function that determines the type of segmentation
% contained in each of the fields. It is used by ft_datatype_segmentation and
% ft_datatype_parcellation.
%
% See also FIXSEGMENTATION, CONVERT_SEGMENTATIONSTYLE

indexed       = false(size(fn));
probabilistic = false(size(fn));

% determine for each of the fields whether it is probabilistic, indexed or something else
for i=1:numel(fn)
  if numel(segmentation.(fn{i}))~=prod(dim)
    % this does not look like a segmentation
    continue
  elseif strcmp(fn{i}, 'anatomy') || strcmp(fn{i}, 'posclusterslabelmat') || strcmp(fn{i}, 'negclusterslabelmat'),
    % these should not be interpreted as segmentation, also not when it is a uint8 or uint16 representation
    continue
  elseif iscell(segmentation.(fn{i}))
    % this should not be interpreted as segmentation
    continue
  else
    if isfield(segmentation, [fn{i} 'label'])
      % the xxxlabel field exists, which only makes sense for an indexed representation
      probabilistic(i) = false;
      indexed(i)       = true;
    else
      % this looks like a segmentation
      tmp = segmentation.(fn{i});
      tmp = tmp(:);       % convert to vector
      sel = isnan(tmp);   % find NaN values
      if any(sel)
        % note that the the finding and removing of NaNs have been separated to speed up the code
        tmp = tmp(~sel);  % remove NaN values
      end
      clear sel
      probabilistic(i) =  islogical(tmp) || all(tmp>=-0.001   & tmp<=1.001); % allow some roundoff error
      indexed(i)       = ~islogical(tmp) && all(tmp>=-0.001) && all(abs(tmp - round(tmp))<1000*eps);
      
      if probabilistic(i) && indexed(i)
        % the xxxlabel does not exist, so treat it as a probabilistic representation
        probabilistic(i) = true;
        indexed(i)       = false;
      end
    end
  end
end % for each of the fields
