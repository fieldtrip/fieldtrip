function segmentation = convert_segmentationstyle(segmentation, fn, dim, style)

% CONVERT_SEGMENTATIONSTYLE is a helper function for converting between probabilistic
% and indexed representations. It is used by ft_datatype_segmentation and
% ft_datatype_parcellation.
%
% See also FIXSEGMENTATION, DETERMINE_SEGMENTATIONSTYLE

switch style
  case 'indexed'
    
    % convert the probabilistic representation to an indexed representation
    threshold = 0.5;
    seg       = zeros(dim);
    seglabel  = cell(1,numel(fn));
    
    for i=1:length(fn)
      tmp = segmentation.(fn{i})>threshold;
      if any(seg(tmp(:)))
        error('overlapping tissue probability maps cannot be converted to an indexed representation');
        % FIXME in principle it is possible to represent two tissue types at one voxel
      else
        seg(tmp) = i;
        seglabel{i} = fn{i};
      end
    end
    segmentation          = rmfield(segmentation, fn);
    segmentation.seg      = seg;
    segmentation.seglabel = seglabel;
    
  case 'probabilistic'
    % convert the indexed representation to a probabilistic one
    for i=1:length(fn)
      fprintf('converting %s\n', fn{i});
      seg      = segmentation.(fn{i});
      seglabel = segmentation.([fn{i} 'label']);
      for j=i:length(seglabel)
        fprintf('creating probabilistic representation for %s\n', seglabel{j});
        segmentation.(fixname(seglabel{j})) = (seg==j);
      end % for j
      segmentation = rmfield(segmentation,  fn{i}         );
      segmentation = rmfield(segmentation, [fn{i} 'label']);
    end % for i
    
  otherwise
    error('unsupported style "%s"', style);
end
