%CPD_AFFINE The affine CPD point-set registration. It is recommended to use
%   and umbrella function rcpd_register with an option opt.method='affine'
%   instead of direct use of the current funciton.
%
%
%   Input
%   ------------------
%   X, Y       real, double, full 2-D matrices of point-set locations. We want to
%              align Y onto X. [N,D]=size(X), where N number of points in X,
%              and D is the dimension of point-sets. Similarly [M,D]=size(Y).
%
%   max_it          maximum number of iterations, try 150
%   tol             tolerance criterium, try 1e-5
%   viz=[0 or 1]    Visualize every iteration         
%   outliers=[0..1] The weight of noise and outliers, try 0.1
%   fgt=[0 or 1]    (default 0) Use a Fast Gauss transform (FGT). (use only for the large data problems)
%   corresp=[0 or 1](default 0) estimate the correspondence vector.
%
%
%
%
%   Output
%   ------------------
%   C      Correspondance vector, such that Y corresponds to X(C,:).
%   B      Affine matrix.
%   t      Translation vector.
%   sigma2 Final sigma^2
%   iter   Final number or iterations
%   T      Registered Y point set
%
%
%   Examples
%   --------
%   It is recommended to use an umbrella function cpd_register with an
%   option opt.method='affine' instead of direct use of the current
%   funciton.
%
%   See also CPD_REGISTER

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

function [C, B, t, sigma2, iter, T]=cpd_affine(X,Y, max_it, tol, viz, outliers, fgt, corresp, sigma2)

[N, D]=size(X);[M, D]=size(Y);

% Initialize sigma and Y
if ~exist('sigma2','var') || isempty(sigma2) || (sigma2==0), 
    sigma2=(M*trace(X'*X)+N*trace(Y'*Y)-2*sum(X)*sum(Y)')/(M*N*D);
end
sigma2_init=sigma2;

T=Y;

% Optimization
iter=0; ntol=tol+10; L=1;
while (iter<max_it) && (ntol > tol) && (sigma2 > 10*eps)

    L_old=L;

    % Check wheather we want to use the Fast Gauss Transform
    if (fgt==0)  % no FGT
        [P1,Pt1, PX, L]=cpd_P(X,T, sigma2 ,outliers); st='';
    else         % FGT
        [P1, Pt1, PX, L, sigma2, st]=cpd_Pfast(X, T, sigma2, outliers, sigma2_init, fgt);
    end
    
    ntol=abs((L-L_old)/L);
    disp([' CPD Affine ' st ' : dL= ' num2str(ntol) ', iter= ' num2str(iter) ' sigma2= ' num2str(sigma2)]);
  
    % Precompute 
    Np=sum(P1);
    mu_x=X'*Pt1/Np;
    mu_y=Y'*P1/Np;


    % Solve for parameters
    B1=PX'*Y-Np*(mu_x*mu_y');
    B2=(Y.*repmat(P1,1,D))'*Y-Np*(mu_y*mu_y');
    B=B1/B2; % B= B1 * inv(B2);
    t=mu_x-B*mu_y;
    
    sigma2save=sigma2;
    sigma2=abs(sum(sum(X.^2.*repmat(Pt1,1,D)))- Np*(mu_x'*mu_x) -trace(B1*B'))/(Np*D); 
    % abs here to prevent roundoff errors that leads to negative sigma^2 in
    % rear cases
    
    % Update centroids positioins
    T=Y*B'+repmat(t',[M 1]);

    iter=iter+1;

    if viz, 
        cpd_plot_iter(X, T); 
%         view(0,-90)
%         axis off
%         fig = gcf;
%         fig.InvertHardcopy = 'off';
%         saveas(fig,['picture',int2str(iter)],'jpg')
    end;
end

% Find the correspondence, such that Y(C) corresponds to X
if corresp, C=cpd_Pcorrespondence(X,T,sigma2save,outliers); else C=0; end;


