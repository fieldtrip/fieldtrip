function [On,Wr] = reorder_mod(W,M)
%REORDER_MOD         Reorder connectivity matrix by modular structure
%
%   On = reorder_mod(W,M);
%   [On Wr] = reorder_mod(W,M);
%
%   This function reorders the connectivity matrix by modular structure and
%   may consequently be useful in visualization of modular structure.
%
%   Inputs:     
%       W,      connectivity matrix (binary/weighted undirected/directed)
%       M,      module affiliation vector
%
%   Outputs:    
%       On,     new node order
%       Wr,     reordered connectivity matrix
%   
%
%   Used in: Rubinov and Sporns (2011) NeuroImage; Zingg et al. (2014) Cell.
%
%
%   2011, Mika Rubinov, UNSW/U Cambridge

%   Modification History:
%   Mar 2011: Original
%	Jan 2015: Improved behavior for directed networks

%#ok<*ASGLU> 
%#ok<*AGROW>

W = W+eps;
[u,dum,M] = unique(M);                                      %make consecutive;
n = numel(M);                                               %number of nodes
m = numel(u);                                               %number of modules

Nm=zeros(1,m);                                              %number of nodes in modules
Knm_o=zeros(n,m);                                           %node-to-module out-degree
Knm_i=zeros(n,m);                                           %node-to-module in-degree
for i=1:m
	Nm(i)=nnz(M==i);
    Knm_o(:,i)=sum(W(:,M==i),2);
    Knm_i(:,i)=sum(W(M==i,:),1);
end
Knm=(Knm_o+Knm_i)/2;

Wm=zeros(m);
for u=1:m
    for v=1:m
        Wm(u,v)=sum(sum(W(M==u,M==v)));
    end
end
Bm=(Wm+Wm.')./(2*(Nm.'*Nm));
% Km_o=sum(Wm,2);
% Km_i=sum(Wm,1);
% Bm_oi=Wm-Km_o*Km_i/sum(sum(Wm));
% Bm=(Bm_oi+Bm_oi.')/2;

%1. Arrange densely connected modules together
[I,J,bv]=find(tril(Bm,-1));                             	%symmetrized intermodular connectivity values
[~,ord]=sort(bv,'descend');                                 %sort by greatest relative connectivity
I=I(ord); 
J=J(ord);
Om=[I(1) J(1)];                                             %new module order

Si=true(size(I));
Sj=true(size(J));
Si(I==I(1) | I==J(1))=0;
Sj(J==I(1) | J==J(1))=0;

while length(Om)<m                                          %while not all modules ordered
    for u=1:numel(I)
        if Si(u)     && any(J(u)==Om([1 end]))
            old=J(u);
            new=I(u);
            break
        elseif Sj(u) && any(I(u)==Om([1 end]))
            old=I(u);
            new=J(u);
            break
        end
    end
    if old==Om(1);      
        Om=[new Om];             
    elseif old==Om(end);
        Om=[Om new];
    end
    Si(I==new)=false;
    Sj(J==new)=false;
end


%2. Reorder nodes within modules
On = zeros(n,1,'uint64');                                   %node order array
for i=1:m
    u = Om(i);
    ind = find(M==u);                                       %indices
        
    mod_imp=[Om; sign((1:m)-i); abs((1:m)-i); Bm(u,Om)].';
    mod_imp=sortrows(mod_imp,[3 -4]);    
    mod_imp=prod(mod_imp(2:end,[1 2]),2);
    
	[dum,ord]=sortrows(Knm(ind,:),mod_imp);                 %sort nodes by number of links to close modules
	On(ind(ord))=1e6*i+(1:Nm(u));                           %assign node order (assumes <1e6 nodes in a module)
end

[dum,On]=sort(On);                                          %reorder nodes
Wr=W(On,On);                                                %reorder matrix