function [R,Nk,Ek] = rich_club_bd(CIJ,varargin)
%RICH_CLUB_BD        Rich club coefficients (binary directed graph)
%
%   R = rich_club_bd(CIJ)
%   [R,Nk,Ek] = rich_club_bd(CIJ,klevel)
%
%   The rich club coefficient, R, at level k is the fraction of edges that
%   connect nodes of degree k or higher out of the maximum number of edges
%   that such nodes might share.
%
%   Input:      CIJ,        connection matrix, binary and directed
%            klevel,        optional input argument. klevel sets the
%                              maximum level at which the rich club
%                              coefficient will be calculated. If klevel is
%                              not included the the maximum level will be
%                              set to the maximum degree of CIJ.
%
%   Output:       R,        vector of rich-club coefficients for levels
%                              1 to klevel.
%                Nk,        number of nodes with degree>k
%                Ek,        number of edges remaining in subgraph with
%                              degree>k
%
%   Reference: Colizza et al. (2006) Nat. Phys. 2:110.
%
%   Martijn van den Heuvel, University Medical Center Utrecht, 2011

N = size(CIJ,1);                    %#ok<NASGU>

% definition of "degree" as used for RC coefficients
% degree is taken to be the sum of incoming and outgoing connectons
[id,od,degree] = degrees_dir(CIJ);  %#ok<ASGLU>

if nargin == 1
    klevel = max(degree);
elseif nargin == 2
    klevel = varargin{1};
elseif nargin > 2
    error('number of inputs incorrect. Should be [CIJ], or [CIJ, klevel]')
end

R = zeros(1,klevel);
Nk = zeros(1,klevel);
Ek = zeros(1,klevel);
for k = 1:klevel
    SmallNodes=find(degree<=k);       %get 'small nodes' with degree <=k
    subCIJ=CIJ;                       %extract subnetwork of nodes >k by removing nodes <=k of CIJ
    subCIJ(SmallNodes,:)=[];          %remove rows
    subCIJ(:,SmallNodes)=[];          %remove columns
    Nk(k)=size(subCIJ,2);             %number of nodes with degree >k
    Ek(k)=sum(subCIJ(:));             %total number of connections in subgraph
    R(k)=Ek(k)/(Nk(k)*(Nk(k)-1));     %unweighted rich-club coefficient
end