function headmodel = ft_headmodel_hbf(mesh, varargin)

% FT_HEADMODEL_HBF creates a volume conduction model of the head
% using the boundary element method (BEM) for M/EEG. This function
% takes as input the triangulated surfaces that describe the boundaries
% and returns as output a volume conduction model which can be used
% to compute leadfields.
% 
% The implementation of this function is based on Matti Stenroos' public
% version of the Hesinki BEM Freamework hence the name "hbf".
%
% Use as
%   headmodel = ft_headmodel_hbf(mesh, ...)
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   conductivity     = 2 x n array, conductivity on inside and outside of
%                               compartment
%   checkmesh        = 'yes' or 'no'
%   isolatedsource   = which compartment number should ISA apply to? 
%                       DEFAULT: []
%
% See also FT_PREPARE_HEADMODEL, FT_COMPUTE_LEADFIELD

% add toolbox
ft_hastoolbox('hbf', 1);

% check if the toolbox has been initialised
if isempty(which('hbf_TM_Phi_LC'))
    hbf_SetPaths;
end
    
% get the optional arguments
conductivity    = ft_getopt(varargin, 'conductivity');
checkmesh       = ft_getopt(varargin, 'checkmesh', 'yes');
isolated        = ft_getopt(varargin, 'isolatedsource', []);

% convert to Boolean value
checkmesh = istrue(checkmesh);

% conductivity values inside and outside of tissue types
ntissues = numel(mesh);
if sum(size(conductivity) == [2 ntissues]) ~= 2
    ft_error('conductivity needs to be an array of size 2-by-%d',ntissues)
end
ci = conductivity(1,:);
co = conductivity(2,:);

% convert meshes into hbf compatible versions
bmeshes = bmesh2bnd(mesh);

% Mesh checking
%--------------

% Perform compulsary integrity check - are the triangles correctly
% oriented? 
for ii = 1:ntissues
    ori = hbf_CheckTriangleOrientation(bmeshes{ii});
    if ori == 2
        fprintf('Changing orientation of triangles in mesh %i\n',ii);
        [bmeshes{ii}, pass] = hbf_CorrectTriangleOrientation(bmeshes{ii});
        if ~pass
            ft_error('Cannot reorient the mesh faces.');
        end
    elseif ori < 1
        ft_error('Cannot determine the orientation of mesh faces.');
    end
end

% Optional mesh checking, on by default
if checkmesh
    % perform general integrity checks
    for ii = 1:ntissues
        status = hbf_CheckMesh(bmeshes{ii});
        if ~status.success
            disp(status)
            msg = ['Mesh %i failed integrity checks. ',...
                'See output above for failed tests'];
            ft_error(msg,ii);
        end

    end

    % if the meshes are nested (for a multi-shell BEM), sort.
    bmeshes = hbf_SortNestedMeshes(bmeshes);
else
    fprintf('Bypassing additional mesh checks\n')
end

% BEM Calculation
%----------------

% Generate transfer matrix
D = hbf_BEMOperatorsPhi_LC(bmeshes);
if isempty(isolated)
    Tphi_full = hbf_TM_Phi_LC(D,ci,co);
else
    fprintf('%-40s: %30s\n','Applying ISA to compartment',...
        num2str(isolated));
    Tphi_full = hbf_TM_Phi_LC_ISA2(D,ci,co,isolated);
end

% Pack up
%--------

headmodel           = [];
headmodel.bnd       = bmesh2bnd(bmeshes); % convert back for plotting
headmodel.cond      = conductivity;
headmodel.mat       = Tphi_full;
headmodel.type      = 'hbf';
headmodel.unit      = mesh(1).unit;

end
