function this = gifti(varargin)
% GIfTI Geometry file format class
% Geometry format under the Neuroimaging Informatics Technology Initiative
% (NIfTI):
%                 http://www.nitrc.org/projects/gifti/
%                      http://nifti.nimh.nih.gov/
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id$

switch nargin
    
    case 0
        this = giftistruct;
        this = class(this,'gifti');
        
    case 1
        if isa(varargin{1},'gifti')
            this = varargin{1};
            
        elseif isstruct(varargin{1})
            f       = {'faces', 'face', 'tri' 'vertices', 'vert', 'pnt', 'cdata'};
            ff      = {'faces', 'faces', 'faces', 'vertices', 'vertices', 'vertices', 'cdata'};
            [c, ia] = intersect(f,fieldnames(varargin{1}));
            if ~isempty(c)
                this = gifti;
                for i=1:length(c)
                    this = subsasgn(this,...
                        struct('type','.','subs',ff{ia(i)}),...
                        varargin{1}.(c{i}));
                end
            elseif isempty(setxor(fieldnames(varargin{1}),...
                    {'metadata','label','data'}))
                this = class(varargin{1},'gifti');
            else
                error('[GIFTI] Invalid structure.');
            end
            
        elseif ishandle(varargin{1})
            this = struct('vertices',get(varargin{1},'Vertices'), ...
                          'faces',   get(varargin{1},'Faces'));
            if ~isempty(get(varargin{1},'FaceVertexCData'));
                  this.cdata = get(varargin{1},'FaceVertexCData');
            end
            this = gifti(this);
            
        elseif isnumeric(varargin{1})
            this = gifti;
            this = subsasgn(this,...
                struct('type','.','subs','cdata'),...
                varargin{1});
            
        elseif ischar(varargin{1})
            if size(varargin{1},1)>1
                this = gifti(cellstr(varargin{1}));
                return;
            end
            [p,n,e] = fileparts(varargin{1});
            if strcmpi(e,'.mat')
                try
                    this = gifti(load(varargin{1}));
                catch
                    error('[GIFTI] Loading of file %s failed.', varargin{1});
                end
            elseif strcmpi(e,'.asc') || strcmpi(e,'.srf')
                this = read_freesurfer_file(varargin{1});
                this = gifti(this);
            else
                this = read_gifti_file(varargin{1},giftistruct);
                this = class(this,'gifti');
            end
            
        elseif iscellstr(varargin{1})
            fnames = varargin{1};
            this(numel(fnames)) = giftistruct;
            this = class(this,'gifti');
            for i=1:numel(fnames)
                this(i) = gifti(fnames{i});
            end
            
        else
            error('[GIFTI] Invalid object construction.');
        end
        
    otherwise
        error('[GIFTI] Invalid object construction.');
end

%==========================================================================
function s = giftistruct
s = struct(...
    'metadata', ...
        struct(...
            'name',       {}, ...
            'value',      {} ...
        ), ...
    'label', ...
        struct(...
            'name',       {}, ...
            'index',      {} ...
        ), ...
    'data', ...
        struct(...
            'attributes', {}, ...
            'metadata',   struct('name',{}, 'value',{}), ...
            'space',      {}, ...
            'data',       {} ...
        ) ...
    );
