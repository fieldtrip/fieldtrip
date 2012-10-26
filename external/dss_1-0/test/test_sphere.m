function test_sphere
%

% Copyright (C) 2004, 2005 DSS MATLAB package team
% (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of
% Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: test_sphere.m,v 1.5 2005/08/29 12:15:52 kosti Exp $

tolerance = 1e-9;

dim = 9;
X = create_mixed_source(dim);
% original source dimension is fixed to 6
sdim = 6;

[W, dW] = dss_sphere(X);

if any(size(W)~=[sdim dim]); error('Invalid sphering matrix dimension'); end
if any(size(dW)~=[dim sdim]); error('Invalid de-sphering matrix dimension'); end
if ~is_identity(W*dW, tolerance); error('Sphering matrices are not inverses'); end
if ~is_identity(cov((W*X)',1), tolerance); error('Sphering matrix doesn''t whiten the data'); end
 
% reduced dimension
rdim = 3;

[W, dW] = dss_sphere(X, rdim);
if any(size(W)~=[rdim dim]); error('Invalid reduced sphering matrix dimension'); end
if any(size(dW)~=[dim rdim]); error('Invalid reduced de-sphering matrix dimension'); end
if ~is_identity(W*dW, tolerance); error('Sphering matrices are not inverses'); end
if ~is_identity(cov((W*X)',1), tolerance); error('Sphering matrix doesn''t whiten the data'); end

[W, dW] = dss_sphere(X,0);
if any(size(W)~=[sdim dim]); error('Dimension parameter 0 reduces dimension'); end
if any(size(dW)~=[dim sdim]); error('Dimension parameter 0 reduces dimension'); end

% try to sphere sphered data
[W2, dW2] = dss_sphere(W*X, rdim);
if ~is_identity(W2,0) | ~is_identity(dW2,0); error('Sphering pre-sphered data returns non-identity sphering matrices.'); end
% pre-sphered data skips dimension reduction
if any(size(W)~=[sdim dim]); error('Invalid pre-sphered reduced sphering matrix dimension'); end
if any(size(dW)~=[dim sdim]); error('Invalid pre-sphered reduced de-sphering matrix dimension'); end

% symmetric whitening
%[W, dW] = sphere(X,0,1);

%if any(size(W)~=[sdim dim]); error('Invalid symmetric sphering matrix dimension'); end
%if any(size(dW)~=[dim sdim]); error('Invalid symmetric de-sphering matrix dimension'); end
%if ~is_identity(W*dW); error('Sphering matrices are not inverses'); end
%if ~is_identity(cov((W*X)',1)); error('Sphering matrix doesn''t whiten the data'); end

function i = is_identity(M, tolerance)
i = ((size(M,1)==size(M,2)) & all(all(abs(M-diag(ones(size(M,1),1))) <= tolerance)));





