%CPD_GRBF_lowrank The nonrigid CPD point-set registration with lowrank
% matrix approximation of G. Works similar to cpd_grbf, but approximate
% Gaussian affinity matrix G by a predefined number of major eigenvectors,
% which allows utilization of Woodburry furmula for fast linear system
% solution.
%   
%   It is recommended to use and umbrella function cpd_register with an option opt.method='nonrigid_lowrank'
%   instead of direct use of the current funciton.
%
%   [C, W, sigma2, iter,T] =cpd_GRBF_lowrank(X, Y, beta, lambda, max_it, tol, viz, outliers, fgt, numeig, corresp);
%
%   Input
%   ------------------
%   X, Y       real, double, full 2-D matrices of point-set locations. We want to
%              align Y onto X. [N,D]=size(X), where N number of points in X,
%              and D is the dimension of point-sets. Similarly [M,D]=size(Y).
%   beta       (try 1..5 ) std of Gaussian filter (Green's funciton) 
%   lambda     (try 1..5) regularization weight
%   max_it     (try 150) maximum number of iterations allowed
%   tol        (try 1e-5) tolerance criterium
%   viz=[0 or 1]    Visualize every iteration         
%   outliers=[0..1] The weight of noise and outliers, try 0.1
%   fgt=[0 or 1]    (default 0) Use a Fast Gauss transform (FGT) for 
%                    a) matrix vectors comutations and b) to find the largest eigenvectors.
%                   (use only for the large data problems)
%   numeig          The number of largest eigenvectors to use, try N^(1/2)
%                   M^(1/2), where M is the length of Y.
%   corresp=[0 or 1](default 0) estimate the correspondence vector.
%
%
%   Output
%   ------------------
%   C      Correspondance vector, such that Y corresponds to X(C,:).
%   W      Non-rigid transformation cooeficients.
%   sigma2 Final sigma^2
%   iter   Final number or iterations
%   T      Registered Y point set
%
%
%   Examples
%   --------
%   It is recommended to use an umbrella function cpd_register with an
%   option opt.method='nonrigid' instead of direct use of the current
%   funciton.
%
%   See also CPD_REGISTER, CPD_GRBF

% Copyright (C) 2008 Andriy Myronenko (myron@csee.ogi.edu)
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


function  [C, W, sigma2, iter, T] =cpd_GRBF_lowrank(X, Y, beta, lambda, max_it, tol, viz, outliers, fgt, numeig, eigfgt, corresp, sigma2);

[N, D]=size(X); [M, D]=size(Y);

T=Y; iter=0;  ntol=tol+10; W=zeros(M,D);
if ~exist('sigma2','var') || isempty(sigma2) || (sigma2==0), 
    sigma2=(M*trace(X'*X)+N*trace(Y'*Y)-2*sum(X)*sum(Y)')/(M*N*D);
end
sigma2_init=sigma2;

% Find 'numeig' eigenvectors of the affinity matrix G
[Q,S]=cpd_GRBF_lowrankQS(Y, beta, numeig, eigfgt);
invS=spdiags(1./diag(S),0,numeig,numeig);


iter=0; ntol=tol+10; L=1;
while (iter<max_it) && (ntol > tol) && (sigma2 > 1e-8)

    
    L_old=L;
    QtW=Q'*W;
    % Check wheather we want to use the Fast Gauss Transform
    if (fgt==0)  % no FGT
        [P1,Pt1, PX, L]=cpd_P(X,T, sigma2 ,outliers); st='';
    else         % FGT
        [P1, Pt1, PX, L, sigma2, st]=cpd_Pfast(X, T, sigma2, outliers, sigma2_init, fgt);
    end
    
    L=L+lambda/2*trace(QtW'*S*QtW);
    ntol=abs((L-L_old)/L);
    disp([' CPD nonrigid ' st ' (lowrank)  : dL= ' num2str(ntol) ', iter= ' num2str(iter) ' sigma2= ' num2str(sigma2)]);

    % M-step. Solve linear system for W.
    dP=spdiags(P1,0,M,M); % precompute diag(P)
  
    dPQ=dP*Q;
    F=PX-dP*Y;

    W=1/(lambda*sigma2)*(F-dPQ*  ( (lambda*sigma2*invS+Q'*dPQ) \ (Q'*F) ) );

    % update Y postions
    % same as T=Y+G*W for full rank;
    T=Y+(Q*(S*(Q'*W)));

    Np=sum(P1);sigma2save=sigma2;
    sigma2=abs(sum(sum(X.^2.*repmat(Pt1,1,D)))+sum(sum(T.^2.*repmat(P1,1,D))) -2*trace(PX'*T) )/(Np*D);
    
    % Plot the result on current iteration
    if viz, cpd_plot_iter(X, T); end;

    iter=iter+1;

end

disp('CPD registration succesfully completed.');

%Find the correspondence, such that Y corresponds to X(C,:)
if corresp, C=cpd_Pcorrespondence(X,T,sigma2save,outliers); else C=0; end;

