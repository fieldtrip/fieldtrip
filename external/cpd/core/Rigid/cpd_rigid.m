%CPD_RIGID The rigid CPD point-set registration. It is recommended to use
%   and umbrella function rcpd_register with an option opt.method='rigid'
%   instead of direct use of the current funciton.
%
%
%   Input
%   ------------------
%   X, Y       real, double, full 2-D matrices of point-set locations. We want to
%              align Y onto X. [N,D]=size(X), where N number of points in X,
%              and D is the dimension of point-sets. Similarly [M,D]=size(Y).
%
%   rot=[0 or 1]    (default 1) 1- strict rotation, 0- allow reflections.
%   scale=[0 or 1]  (default 1) 1- estimate scaling, 0- don't estimate scaling.    
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
%   R      Rotation matrix.
%   t      Translation vector.
%   s      Scaling constant.
%   sigma2 Final sigma^2
%   iter   Final number or iterations
%   T      Registered Y point set
%
%
%   Examples
%   --------
%   It is recommended to use an umbrella function cpd_register with an
%   option opt.method='rigid' instead of direct use of the current
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

function [C, R, t, s, sigma2, iter, T]=cpd_rigid(X,Y, rot, scale, max_it, tol, viz, outliers, fgt, corresp, sigma2)

[N, D]=size(X);[M, D]=size(Y);
if viz, figure; end;

% Initialization
if ~exist('sigma2','var') || isempty(sigma2) || (sigma2==0), 
    sigma2=(M*trace(X'*X)+N*trace(Y'*Y)-2*sum(X)*sum(Y)')/(M*N*D);
end
sigma2_init=sigma2;

% % add prior information
% reliability = 1;
% p_prior     = zeros(length(prior1),length(prior2));
% for i = length(prior1)
%     for j = length(prior2)
%         p_prior(i,j) = 1/(2*pi*reliability)^(3/2)*exp(-norm(prior1(i,1:3)-prior2(j,1:3),2)/(2*reliability^2));
%     end    
% end


T=Y; s=1; R=eye(D);
% Optimization
iter=0; ntol=tol+10; L=0; Y_old = Y;checker = 1;
while (iter<max_it) && (ntol > tol) && (sigma2 > 10*eps)

    L_old=L;
    % Check wheather we want to use the Fast Gauss Transform
    if (fgt==0)  % no FGT
        [P1,Pt1, PX, L]=cpd_P(X,T, sigma2 ,outliers); st='';
    else         % FGT
        [P1, Pt1, PX, L, sigma2, st]=cpd_Pfast(X, T, sigma2, outliers, sigma2_init, fgt);
    end

    ntol=abs((L-L_old)/L);
    %disp([' CPD Rigid ' st ' : dL= ' num2str(ntol) ', iter= ' num2str(iter) ' sigma2= ' num2str(sigma2)]);

    % Precompute
    Np=sum(Pt1);
    mu_x=X'*Pt1/Np;
    mu_y=Y'*P1/Np;

    % Solve for Rotation, scaling, translation and sigma^2
    A=PX'*Y-Np*(mu_x*mu_y'); % A= X'P'*Y-X'P'1*1'P'Y/Np;
    [U,S,V]=svd(A); C=eye(D);
    if rot, C(end,end)=det(U*V'); end % check if we need strictly rotation (no reflections)
    R=U*C*V';

    sigma2save=sigma2;
    if scale  % check if estimating scaling as well, otherwise s=1
        s=trace(S*C)/(sum(sum(Y.^2.*repmat(P1,1,D))) - Np*(mu_y'*mu_y));
        sigma2=abs(sum(sum(X.^2.*repmat(Pt1,1,D))) - Np*(mu_x'*mu_x) -s*trace(S*C))/(Np*D);
    else
        sigma2=abs((sum(sum(X.^2.*repmat(Pt1,1,D))) - Np*(mu_x'*mu_x)+sum(sum(Y.^2.*repmat(P1,1,D))) - Np*(mu_y'*mu_y) -2*trace(S*C))/(Np*D));
    end

    t=mu_x-s*R*mu_y;

    % Update the GMM centroids
    T=s*Y_old*R'+repmat(t',[M 1]);

    iter=iter+1;

    if viz
        cpd_plot_iter(X, T);  
%         view([0 -90]);
%         xlim([-1.5 1.5])
%         ylim([-1.5 1.5])
%         zlim([-1.5 1.5])
%         axis off
%         filename = sprintf('test_image_%03d.jpeg', iter);
%         saveas(gcf, filename, 'jpeg') ;
    end % show current iteration if viz=1
end
% Find the correspondence, such that Y corresponds to X(C,:)
if corresp, C=cpd_Pcorrespondence(X,T,sigma2save,outliers); else C=0; end;



