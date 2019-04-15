function this = subsasgn(this, subs, A)
% Subscript assignment for GIfTI objects
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: subsasgn.m 6513 2015-08-05 17:52:13Z guillaume $

switch subs(1).type
    case '.'
        if ~ismember(subs(1).subs, ...
                {'vertices' 'faces' 'normals' 'cdata','mat','indices','private'})
            error('Reference to non-existent field ''%s''.',subs(1).subs);
        end
        % TODO % handle cases when length(subs) > 1
        [i,n] = isintent(this,subs(1).subs);
        if isempty(i) && ~strcmp(subs(1).subs,'private')
            n = length(this.data) + 1;
            if n==1, this.data = {}; end
            % TODO % Initialise data field appropriately
            this.data{n}.metadata = struct([]);
            this.data{n}.space    = [];
            this.data{n}.attributes.Dim = size(A);
            % TODO % set DataType according to intent type
            this.data{n}.data = [];
            switch subs(1).subs
                case {'vertices','mat'}
                    in = 'NIFTI_INTENT_POINTSET';
                    dt = 'NIFTI_TYPE_FLOAT32';
                    this.data{n}.space.DataSpace = 'NIFTI_XFORM_UNKNOWN';
                    this.data{n}.space.TransformedSpace = 'NIFTI_XFORM_UNKNOWN';
                    this.data{n}.space.MatrixData = eye(4);
                case 'faces'
                    in = 'NIFTI_INTENT_TRIANGLE';
                    dt = 'NIFTI_TYPE_INT32';
                case 'indices'
                    in = 'NIFTI_INTENT_NODE_INDEX';
                    dt = 'NIFTI_TYPE_INT32';
                case 'normals'
                    in = 'NIFTI_INTENT_VECTOR';
                    dt = 'NIFTI_TYPE_FLOAT32';
                case 'cdata'
                    in = 'NIFTI_INTENT_NONE';
                    dt = 'NIFTI_TYPE_FLOAT32';
                otherwise
                    error('This should not happen.');
            end
            this.data{n}.attributes.Intent = in;
            this.data{n}.attributes.DataType = dt;
        end
        
        switch subs(1).subs            
            %- .private
            %--------------------------------------------------------------
            case 'private'
                this = builtin('subsasgn',this,subs(2:end),A);
                
            % .mat
            %--------------------------------------------------------------
            case 'mat'
                if length(subs) > 1
                    this.data{n}.space(1).MatrixData = builtin('subsasgn',...
                        this.data{n}.space(1).MatrixData,subs(2:end),A);
                else
                    if ~isequal(size(A),[4 4])
                        error('Invalid Coordinate System Transform Matrix.');
                    end
                    this.data{n}.space(1).MatrixData = A;
                end
                
            %- .faces
            %--------------------------------------------------------------
            case 'faces'
                if length(subs) > 1
                    this.data{n}.data = int32(builtin('subsasgn',this.data{n}.data,subs(2:end),A-1));
                else
                    this.data{n}.data = int32(A - 1);
                    this.data{n}.attributes.Dim = size(A);
                end
                
            %- .indices
            %--------------------------------------------------------------
            case 'indices'
                if n ~= 1
                    this.data = this.data([n setdiff(1:numel(this.data),n)]);
                    n = 1;
                end
                if length(subs) > 1
                    this.data{n}.data = int32(builtin('subsasgn',this.data{n}.data,subs(2:end),A-1));
                else
                    A = A(:);
                    this.data{n}.data = int32(A - 1);
                    this.data{n}.attributes.Dim = size(A);
                end
                
            %- .vertices, .normals, .cdata
            %--------------------------------------------------------------
            otherwise
                if length(subs) > 1
                    if numel(n) == 1
                        this.data{n}.data = single(builtin('subsasgn',this.data{n}.data,subs(2:end),A));
                        this.data{n}.attributes.Dim = size(this.data{n}.data);
                    else
                        if numel(subs(2).subs) == 1
                            error('Linear indexing not supported: use multiple subscripts.');
                        end
                        idx = subs(2).subs{2};
                        if isequal(idx,':'), idx = 1:numel(this.data); end
                        for k=1:numel(idx)
                            s = subs(2);
                            s.subs{2} = 1;
                            if numel(A) == 1
                                this.data{idx(k)}.data = single(builtin('subsasgn',this.data{idx(k)}.data,s,A));
                            else
                                this.data{idx(k)}.data = single(builtin('subsasgn',this.data{idx(k)}.data,s,A(:,k)));
                            end
                            this.data{idx(k)}.attributes.Dim = size(this.data{idx(k)}.data);
                        end
                    end
                else
                    if numel(n) == 1
                        if isa(A,'file_array')
                            this.data{n}.data = A;
                            this.data{n}.attributes.Dim = A.dim;
                        else
                            this.data{n}.data = single(A);
                            this.data{n}.attributes.Dim = size(A);
                        end
                    else
                        error('Syntax not implemented.');
                    end
                end
        end
        
    case '()'
    case '{}'
    otherwise
        error('This should not happen.');
end
