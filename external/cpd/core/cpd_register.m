%CPD_REGISTER Rigid, affine, non-rigid point set registration. 
% The main CPD registration function that sets all the options,
% normalizes the data, runs the rigid/non-rigid registration, and returns
% the transformation parameters, registered poin-set, and the
% correspondences.
%
%   Input
%   ------------------ 
%   X, Y       real, double, full 2-D matrices of point-set locations. We want to
%              align Y onto X. [N,D]=size(X), where N number of points in X,
%              and D is the dimension of point-sets. Similarly [M,D]=size(Y). 
%   
%   opt        a structure of options with the following optional fields:
%
%       .method=['rigid','affine','nonrigid','nonrigid_lowrank'] (default
%               rigid) Registration method. Nonrigid_lowrank uses low rank
%               matrix approximation (use for the large data problems).
%       .corresp=[0 or 1] (default 0) estimate the correspondence vector at
%               the end of the registration.
%       .normalize =[0 or 1] (default 1) - normalize both point-sets to zero mean and unit
%               variance before registration, and denormalize after.
%       .max_it (default 150) Maximum number of iterations.
%       .tol (default 1e-5) Tolerance stopping criterion.
%       .viz=[0 or 1] (default 1) Visualize every every iteration.
%       .outliers=[0...1] (default 0.1) The weight of noise and outliers
%       .fgt=[0,1 or 2] Default 0 - do not use. 1 - Use a Fast Gauss transform (FGT) to speed up some matrix-vector product.
%                       2 - Use FGT and fine tune (at the end) using truncated kernel approximations.
%
%       Rigid registration options
%       .rot=[0 or 1] (default 1) 1 - estimate strictly rotation. 0 - also allow for reflections.
%       .scale=[0 or 1] (default 1) 1- estimate scaling. 0 - fixed scaling. 
%
%       Non-rigid registration options
%       .beta [>0] (default 2) Gaussian smoothing filter size. Forces rigidity.
%       .lambda [>0] (default 3) Regularization weight. 
%
%   Output
%   ------------------ 
%   Transform      structure of the estimated transformation parameters:
%
%           .Y     registered Y point-set
%           .iter  total number of iterations
%
%                  Rigid/affine cases only:     
%           .R     Rotation/affine matrix.
%           .t     Translation vector.
%           .s     Scaling constant.
%           
%                  Non-rigid cases only:
%           .W     Non-rigid coefficient
%           .beta  Gaussian width
%           .t, .s translation and scaling
%           
%    
%   C       Correspondance vector, such that Y corresponds to X(C,:)
%
%           
%
%   Examples
%   --------
%
%   See many detailed examples in the 'examples' folder.

% Copyright (C) 2008-2009 Andriy Myronenko (myron@csee.ogi.edu)
% also see http://www.bme.ogi.edu/~myron/matlab/cpd/
%
%     This file is part of the Coherent Point Drift (CPD) package.
%
%     The source code is provided under the terms of the GNU General Public License as published by
%     the Free Software Foundation version 2 of the License.
% 
%     CPD package is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with CPD package; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

function [Transform, C]=cpd_register(X, Y, opt)


[M,D]=size(Y); [N, D2]=size(X);
% Check the input options and set the defaults
if nargin<2, error('cpd_rigid_register.m error! Not enough input parameters.'); end;
if nargin<3, opt.method='rigid'; end;
if ~isfield(opt,'method') || isempty(opt.method), opt.method = 'rigid'; end;
if ~isfield(opt,'normalize') || isempty(opt.normalize), opt.normalize = 1; end;
if ~isfield(opt,'max_it') || isempty(opt.max_it), opt.max_it = 150; end;
if ~isfield(opt,'tol') || isempty(opt.tol), opt.tol = 1e-5; end;
if ~isfield(opt,'viz') || isempty(opt.viz), opt.viz = 1; end;
if ~isfield(opt,'corresp') || isempty(opt.corresp), opt.corresp = 0; end;
if ~isfield(opt,'outliers') || isempty(opt.outliers), opt.outliers = 0.1; end;
if ~isfield(opt,'fgt') || isempty(opt.fgt), opt.fgt = 0; end;
if ~isfield(opt,'sigma2') || isempty(opt.sigma2), opt.sigma2 = 0; end;

% strictly rigid params
if ~isfield(opt,'rot') || isempty(opt.rot), opt.rot = 1; end;
if ~isfield(opt,'scale') || isempty(opt.scale), opt.scale = 1; end;
% strictly non-rigid params
if ~isfield(opt,'beta') || isempty(opt.beta), opt.beta = 2; end;
if ~isfield(opt,'lambda') || isempty(opt.lambda), opt.lambda = 3; end;
% lowrank app param
if ~isfield(opt,'numeig') || isempty(opt.numeig), opt.numeig = round(sqrt(M)); end;
if ~isfield(opt,'eigfgt') || isempty(opt.eigfgt), opt.eigfgt = 1; end;


% checking for the possible errors
if D~=D2, error('The dimension of point-sets is not the same.'); end;
if (D>M)||(D>N), disp('The dimensionality is larger than the number of points. Possibly the wrong orientation of X and Y.'); end;
if (M>1e+5)||(N>1e+5) && (opt.fgt==0), disp('The data sets are large. Use opt.fgt=1 to speed up the process.'); end;
if (M>1e+5)||(N>1e+5) && strcmp(lower(opt.method),'nonrigid'), disp('The data sets are large. Use opt.method=''nonrigid_lowrank'' to speed up the non-rigid registration'); end;
if (D<=1) || (D>3), opt.viz=0; end;
if (opt.normalize==1)&&(opt.scale==0),opt.scale=1; end;

% check if mex functions are compiled yet
if ~exist('cpd_P','file')
    disp('Looks like you have not compiled CPD mex files yet (needs to be done once)');
    disp('Running cpd_make.m for you ...'); tic;
    cpd_make;
end

% Convert to double type, save Y
X=double(X);  
Y=double(Y); Yorig=Y; 

% default mean and scaling
normal.xd=0; normal.yd=0;
normal.xscale=1; normal.yscale=1;

% Normalize to zero mean and unit variance
if opt.normalize, [X,Y,normal]=cpd_normalize(X,Y); end;

disp(['%%%%% Starting CPD-' upper(opt.method) ' registration. %%%' ]); tic;

%%%% Choose the method and start CPD point-set registration
switch lower(opt.method),
    case 'rigid'
        [C, R, t, s, sigma2, iter, T]=cpd_rigid(X,Y, opt.rot, opt.scale, opt.max_it, opt.tol, opt.viz, opt.outliers, opt.fgt, opt.corresp, opt.sigma2);
       case 'affine'
        [C, R, t, sigma2, iter, T]=cpd_affine(X,Y, opt.max_it, opt.tol, opt.viz, opt.outliers, opt.fgt, opt.corresp, opt.sigma2); s=1;
    case 'nonrigid'
        [C, W, sigma2, iter, T] =cpd_GRBF(X, Y, opt.beta, opt.lambda, opt.max_it, opt.tol, opt.viz, opt.outliers, opt.fgt, opt.corresp, opt.sigma2);    
    case 'nonrigid_lowrank'
        [C, W, sigma2, iter, T] =cpd_GRBF_lowrank(X, Y, opt.beta, opt.lambda, opt.max_it, opt.tol, opt.viz, opt.outliers, opt.fgt, opt.numeig, opt.eigfgt, opt.corresp, opt.sigma2);
    otherwise
        error('The opt.method value is invalid. Supported methods are: rigid, affine, nonrigid, nonrigid_lowrank');
end
%%%% 
disptime(toc);

Transform.iter=iter;
Transform.method=opt.method;
Transform.Y=T;
Transform.normal=normal;

% Denormalize transformation parameters
switch lower(opt.method)
    case {'rigid','affine'}
        Transform.R=R; Transform.t=t;Transform.s=s;
        if opt.normalize, % denormalize parameters and registered point set, if it was prenormalized
            Transform.s=Transform.s*(normal.xscale/normal.yscale);
            Transform.t=normal.xscale*Transform.t+normal.xd'-Transform.s*(Transform.R*normal.yd');
            Transform.Y=T*normal.xscale+repmat(normal.xd,M,1);
            
            if strcmp(lower(opt.method),'affine')
                Transform.R=Transform.s*Transform.R; 
                Transform.s=1;
            end
            
        end;
    case {'nonrigid','nonrigid_lowrank'}
            Transform.beta=opt.beta;
            Transform.W=W;
            Transform.Yorig=Yorig;
            Transform.s=1;
            Transform.t=zeros(D,1);
        if opt.normalize,
             Transform.Y=T*normal.xscale+repmat(normal.xd,M,1);
        end 

end




