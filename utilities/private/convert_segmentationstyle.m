function segmentation = convert_segmentationstyle(segmentation, fn, dim, style)

% CONVERT_SEGMENTATIONSTYLE is a helper function for converting between probabilistic
% and indexed representations. It is used by FT_DATATYPE_SEGMENTATION and
% FT_DATATYPE_PARCELLATION.
%
% See also FIXSEGMENTATION, DETERMINE_SEGMENTATIONSTYLE

switch style
  case 'indexed'
    
    % convert the probabilistic representation to an indexed representation
    threshold   = 0.5;
    tissue      = zeros(dim);
    tissuelabel = cell(1,numel(fn));
    
    for i=1:length(fn)
      tmp = segmentation.(fn{i})>threshold;
      if any(tissue(tmp(:)))
        ft_error('overlapping tissue probability maps cannot be converted to an indexed representation');
        % FIXME in principle it is possible to represent two tissue types at one voxel
      else
        tissue(tmp) = i;
        tissuelabel{i} = fn{i};
      end
    end
    segmentation             = rmfield(segmentation, fn);
    segmentation.tissue      = tissue;
    segmentation.tissuelabel = tissuelabel;
    
  case 'probabilistic'
    % convert the indexed representation to a probabilistic one
    for i=1:length(fn)
      fprintf('converting %s\n', fn{i});
      tissue      = segmentation.(fn{i});
      tissuelabel = segmentation.([fn{i} 'label']);
      for j=i:length(tissuelabel)
        fprintf('creating probabilistic representation for %s\n', tissuelabel{j});
        tmp.(fixname(tissuelabel{j})) = (tissue==j);  % avoid overwriting the existing segmentation
      end % for j
      segmentation = rmfield(segmentation,  fn{i}         );
      segmentation = rmfield(segmentation, [fn{i} 'label']);
      for j=i:length(tissuelabel)
        segmentation.(fixname(tissuelabel{j})) = tmp.(fixname(tissuelabel{j})); % avoid overwriting the existing segmentation
      end
    end % for i
    
  otherwise
    ft_error('unsupported style "%s"', style);
end
