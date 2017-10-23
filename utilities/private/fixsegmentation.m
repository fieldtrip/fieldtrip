function segmentation = fixsegmentation(segmentation, fn, style)

% FIXSEGMENTATION is a helper function that ensures the segmentation to be internally
% consistent. It is used by ft_datatype_segmentation and ft_datatype_parcellation.
%
% % See also CONVERT_SEGMENTATIONSTYLE, DETERMINE_SEGMENTATIONSTYLE

switch style
  case 'indexed'
    
    for i=1:length(fn)
      indexval = unique(segmentation.(fn{i})(:));  % find the unique tissue types
      indexval = indexval(indexval~=0);            % these are the only ones that matter
      
      if any(indexval<0)
        ft_error('an indexed representation cannot contain negative numbers');
      end
      
      if ~isfield(segmentation, [fn{i} 'label'])
        % ensure that the tissues have labels
        indexlabel = cell(size(indexval));
        for j=1:length(indexval)
          indexlabel{indexval(j)} = sprintf('tissue %d', indexval(j));
        end
        segmentation.([fn{i} 'label']) = indexlabel;
        
      else
        % check for the situation where
        %   indexval   = [1 2 4]
        %   indexlabel = {'a', 'b', 'c', 'd'} or {'a', 'b', [], 'd'}
        % which happens if the segmentation unexpectedly does not contain a certain tissue type
        indexlabel = segmentation.([fn{i} 'label']);
        if numel(indexval)>numel(indexlabel)
          ft_error('each index value should have a corresponding entry in %s', [fn{i} 'label']);
        elseif any(cellfun(@isempty, indexlabel(indexval)))
          ft_error('each index value should have a corresponding entry in %s', [fn{i} 'label']);
        elseif numel(indexval)<numel(indexlabel)
          ft_warning('there are more labels than actuall tissue types in %s', [fn{i} 'label']);
          % ensure that the indices are subsequent integers, i.e. [1 2 3] rather than [1 2 4]
          for j=1:length(indexval)
            tmp = segmentation.(fn{i});
            tmp(tmp==indexval(j)) = j;
            segmentation.(fn{i}) = tmp;
          end
          segmentation.([fn{i} 'label']) = segmentation.([fn{i} 'label'])(indexval);
        end
        
      end
    end
    clear tmp indexval indexlabel
    
  case 'probabilistic'
    
    % convert from a cumulative to an exclusive representation
    contains = false(length(fn));
    if length(fn)>4
      % test for each tissue whether it is overlapping with or contained in each other tissue
      ft_warning('more than 4 tissue types, this may take a while');
    end
    for i=1:length(fn)
      segi = segmentation.(fn{i})>0;
      for j=1:length(fn)
        if i==j
          % don't test for self-overlap
          continue
        end
        if ~any(segi(:))
          % don't bother to test completely empty segmentations
          continue
        end
        segj = segmentation.(fn{j})>0;
        contains(i,j) = all(segj(segi(:))); % segi is fully contained in segj
        if i~=j && contains(i,j)
          fprintf('the %s is fully contained in the %s, removing it from the %s\n', fn{i}, fn{j}, fn{j});
          segmentation.(fn{j})(segi) = 0;
        end
      end
    end
    clear segi segj contains
    
  otherwise
    ft_error('unsupported style "%s"', style);
end
