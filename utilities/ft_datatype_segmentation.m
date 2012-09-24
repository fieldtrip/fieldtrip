function segmentation = ft_datatype_segmentation(segmentation, varargin)

% FT_DATATYPE_SEGMENTATION describes the FieldTrip MATLAB structure for
% segmented volumetric data and atlases.
%
% A segmentation is a volumetric description which is usually derived from
% an anatomical MRI, which describes for each voxel the tissue type. It for
% example distinguishes between white matter, grey matter, csf, skull and
% skin. It is mainly used for masking in visualization, construction of
% volume conduction models and for construction of cortical sheets. An
% volume-based atlas is basically a very detailled segmentation with
% anatomical labels for each tissue type.
%
% For example, the AFNI TTatlas+tlrc segmented brain atlas (which can be
% created with FT_PREPARE_ATLAS) looks like this
%
%              dim: [161 191 141]        the size of the 3D volume in voxels
%        transform: [4x4 double]         affine transformation matrix for mapping the voxel coordinates to head coordinate system
%            coord: 'tal'                the transformation matrix maps the voxels into this (head) coordinate system
%             unit: 'mm'                 the units in which the coordinate system is expressed
%           brick0: [161x191x141 uint8]  integer values from 1 to N, the value 0 means unknown
%           brick1: [161x191x141 uint8]  integer values from 1 to M, the value 0 means unknown
%      brick0label: {Nx1 cell}
%      brick1label: {Mx1 cell}
%
% An example of a whole-brain anatomical MRI that was segmented using FT_VOLUMESEGMENT looks like this
%
%         dim: [256 256 256]         the size of the 3D volume in voxels
%   transform: [4x4 double]          affine transformation matrix for mapping the voxel coordinates to head coordinate system
%    coordsys: 'ctf'                 the transformation matrix maps the voxels into this (head) coordinate system
%        unit: 'mm'                  the units in which the coordinate system is expressed
%        gray: [256x256x256 double]  probabilistic map of the gray matter
%       white: [256x256x256 double]  probabilistic map of the white matter
%         csf: [256x256x256 double]  probabilistic map of the cerebrospinal fluid
%
% An example segmentation with binary values that can be used for
% construction of a BEM volume conduction model of the head looks like this
%
%           dim: [256 256 256]         the dimensionality of the 3D volume
%     transform: [4x4 double]          affine transformation matrix for mapping the voxel coordinates to head coordinate system
%      coordsys: 'ctf'                 the transformation matrix maps the voxels into this (head) coordinate system
%          unit: 'mm'                  the units in which the coordinate system is expressed
%         brain: [256x256x256 logical] binary map representing the voxels which belong to the brain
%         scalp: [256x256x256 logical] binary map representing the voxels which belong to the scalp
%         skull: [256x256x256 logical] binary map representing the voxels which belong to the skull
%
% The only difference to the volume data representation is that the segmentation
% structure contains the additional fields XXX. See FT_DATATYPE_VOLUME for
% further details.
%

% get the optional input arguments, which should be specified as key-value pairs
version           = ft_getopt(varargin, 'version', 'latest');
segmentationstyle = ft_getopt(varargin, 'segmentationstyle');   % can be indexed or probabilistic
hasbrainmask      = ft_getopt(varargin, 'hasbrainmask', 'no');  % no means that it is not required, if present it won't be removed

% convert from string into boolean
hasbrainmask = istrue(hasbrainmask);

if strcmp(version, 'latest')
  segversion = '2012';
  volversion = 'latest';
else
  segversion = version;
  volversion = version;
end

if isempty(segmentation)
  return;
end

switch segversion
  case '2012'
    % determine whether the style of the input fields is probabilistic or indexed
    fn = fieldnames(segmentation);
    indexed       = false(size(fn));
    probabilistic = false(size(fn));
    
    [dum, i] = intersect(fn, {'scalp', 'skull', 'brain'});
    if numel(i)==3
      % put them in the preferred order
      fn(i) = {'scalp', 'skull', 'brain'};
    end
    [dum, i] = intersect(fn, {'skin', 'skull', 'brain'}); % this is not likely
    if numel(i)==3
      % put them in the preferred order
      fn(i) = {'skin', 'skull', 'brain'};
    end
    
    % determine for each of the fields whether it is probabilistic, indexed or something else
    for i=1:numel(fn)
      if numel(segmentation.(fn{i}))~=prod(segmentation.dim)
        % this does not look like a segmentation
        continue
      else
        if isfield(segmentation, [fn{i} 'label'])
          % the xxxlabel field exists, which only makes sense for an indexed representation
          probabilistic(i) = false;
          indexed(i)       = true;
        else
          % this looks like a segmentation
          tmp = segmentation.(fn{i});
          tmp = tmp(:);
          probabilistic(i) =  islogical(tmp) || all(tmp>=-0.001 & tmp<=1.001); % allow some roundoff error
          indexed(i)       = ~islogical(tmp) && all(abs(tmp - round(tmp))<1000*eps);
          
          if probabilistic(i) && indexed(i)
            % the xxxlabel does not exist, so treat it as a probabilistic representation
            probabilistic(i) = true;
            indexed(i)       = false;
          end
        end
      end
    end % for each of the fields
    
    if any(probabilistic) && any(indexed)
      warning('cannot work with a mixed representation, removing tissue probability maps');
      sel = find(probabilistic);
      segmentation       = rmfield(segmentation, fn(sel));
      probabilistic(sel) = false;
    end
    
    if any(probabilistic)
      fn = fn(probabilistic);
      probabilistic = true; indexed = false;
    elseif any(indexed)
      fn = fn(indexed);
      indexed = true; probabilistic = false;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fix the indexed representation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if indexed
      for i=1:length(fn)
        indexval = unique(segmentation.(fn{i})(:));  % find the unique tissue types
        indexval = indexval(indexval>0);             % these are the only ones that matter
        % ensure that the tissues have labels
        if ~isfield(segmentation, [fn{i} 'label'])
          indexlabel = {};
          for j=1:length(indexval)
            indexlabel{indexval(j)} = sprintf('tissue%d', indexval(j));
          end
          segmentation.([fn{i} 'label']) = indexlabel;
        end
        % ensure that the indices are subsequent integers, i.e. [1 2 3] rather than [1 2 4]
        for j=1:length(indexval)
          tmp = segmentation.(fn{i});
          tmp(tmp==indexval(j)) = j;
          segmentation.(fn{i}) = tmp;
        end
        segmentation.([fn{i} 'label']) = segmentation.([fn{i} 'label'])(indexval);
      end
      clear tmp indexval indexlabel
    end % if indexed
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fix the probabilistic representation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if probabilistic
      % convert from a cumulative to an exclusive representation
      contains = false(length(fn));
      if length(fn)>4
        % test for each tissue whether it is overlapping with or contained in each other tissue
        warning('more than 4 tissue types, this may take a while');
      end
      for i=1:length(fn)
        segi = segmentation.(fn{i})>0;
        for j=1:length(fn)
          if i==j
            % don't test for self-overlap
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
    end % if probabilistic
    
    % convert from an exclusive to cumulative representation
    % this is only only for demonstration purposes
    % for i=1:length(sel)
    %   segmentation.(fn{sel(i)}) = volumefillholes(segmentation.(fn{sel(i)}));
    % end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert the segmentation to the desired style
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if probabilistic && strcmp(segmentationstyle, 'indexed')
      
      % convert the probabilistic representation to an indexed representation
      threshold = 0.5;
      seg       = zeros(segmentation.dim);
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
      
      clear fn % to avoid any potential confusion
      indexed       = true;
      probabilistic = false;
      
    elseif indexed && strcmp(segmentationstyle, 'probabilistic')
      
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
      
      clear fn % to avoid any potential confusion
      probabilistic = true;
      indexed       = false;
      
    end % converting between probabilistic and indexed
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add the brainmask if requested
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if hasbrainmask
      if indexed
        fn = fieldnames(segmentation);
        sel = false(size(fn));
        for i=1:numel(fn)
          sel(i) = any(strcmp(fn, [fn{i} 'label']));
        end
        fn = fn(sel);
        
        if numel(fn)>1
          error('cannot construct a brain mask on the fly; this requires a single indexed representation');
        else
          seg      = segmentation.(fn{1});
          seglabel = segmentation.([fn{1} 'label']);
          if ~any(strcmp(seglabel, 'brain'))
            threshold = 0.5;
            smooth    = 5;
            % ensure that the segmentation contains the brain mask, if not then construct it from gray+white+csf
            if length(intersect(seglabel, {'gray' 'white' 'csf'}))~=3
              error('cannot construct a brain mask on the fly; this requires gray, white and csf');
            end
            gray  = seg==find(strcmp(seglabel, 'gray'));
            white = seg==find(strcmp(seglabel, 'white'));
            csf   = seg==find(strcmp(seglabel, 'csf'));
            brain = gray + white + csf;
            clear gray white csf seg
            brain = volumesmooth(brain,    smooth,    'brainmask');
            brain = volumethreshold(brain, threshold, 'brainmask');
            % store it in the output
            segmentation.brain = brain;
          end % try to construct the brainmask
        end
        
        
      elseif probabilistic
        if ~isfield(segmentation, 'brain')
          if ~all(isfield(segmentation, {'gray' 'white' 'csf'}))
            error('cannot construct a brain mask on the fly; this requires gray, white and csf');
          end
          threshold = 0.5;
          smooth    = 5;
          % ensure that the segmentation contains the brain mask, if not then construct it from gray+white+csf tissue probability maps
          gray  = segmentation.gray;
          white = segmentation.white;
          csf   = segmentation.csf;
          brain = gray + white + csf;
          clear gray white csf
          brain = volumesmooth(brain,    smooth,    'brainmask');
          brain = volumethreshold(brain, threshold, 'brainmask');
          % store it in the output
          segmentation.brain = brain;
        end
      end
    end % if hasbrainmask
    
  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('unsupported version "%s" for segmentation datatype', version);
end

% convert it into a volume structure
volume = segmentation;

% remove the fields that are specific for the segmentation
% FIXME

volume = ft_datatype_volume(volume, 'version', volversion);

% add the segmentation specific fields again
% FIXME

segmentation     = volume;
clear volume

