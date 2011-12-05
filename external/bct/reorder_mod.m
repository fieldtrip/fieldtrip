function [On Ar] = reorder_mod(A,M)
%REORDER_MOD         Reorder connectivity matrix by modular structure
%
%   On = reorder_mod(A,M);
%
%   This function reorders the connectivity matrix by modular structure and
%   may hence be useful in visualization of modular structure.
%
%   Inputs:     A,          connectivity matrix (with no negative weights)
%               M,          module affiliation vector
%
%   Outputs:    On,         new node order
%               Ar,         reordered connectivity matrix
%   
%
%   Used in: Rubinov and Sporns (2011) NeuroImage.
%
%
%   2011, Mika Rubinov, UNSW

%   Modification History:
%   Mar 2011: Original


[u dum M] = unique(M);                                      %make consecutive;
n = numel(M);                                               %number of nodes
m = numel(u);                                               %number of modules

Nm=zeros(1,m);                                              %number of nodes in modules
Knm=zeros(n,m);                                             %degree to other modules
for i=1:m
	Nm(i)=nnz(M==i);
	Knm(:,i)=sum(A(:,M==i),2);
end
% Kn=sum(A,2);

Am=zeros(m);                                                %relative intermodular connectivity
for i=1:m;
    Am(i,:)=sum(Knm(M==i,:),1);
end
Am=Am./(Nm.'*Nm);

%1. Arrange densely connected modules together
[I J s]=find(tril(Am,-1)+1);                             	%symmetrized intermodular connectivity values
[dum ord]=sort(s,'descend');                                %sort by greatest relative connectivity
I=I(ord); J=J(ord);
Om=[I(1) J(1)];                                             %new module order
% Omq=Om;                                                     %module queue for node reordering
I(1)=0;   J(1)=0;
while length(Om)<m                                          %while not all modules ordered
    ui = find(I & ((J==Om(1)) | J==Om(end)), 1);
    uj = find(J & ((I==Om(1)) | I==Om(end)), 1);
    if ui==uj;  I(ui)=0; J(uj)=0; continue;     end
    
    if isempty(ui);     ui=inf;                 end
    if isempty(uj);     uj=inf;                 end
    if ui<uj;           old=J(ui); new=I(ui);   end
    if uj<ui;           old=I(uj); new=J(uj);   end
    if old==Om(1);      Om=[new Om];            end
    if old==Om(end);    Om=[Om new];            end
    
    I(I==old)=0; J(J==old)=0;
%     Omq=[Omq new];
end

% Nm=Nm(Om);
% Knm=Knm(:,Om);

%2. Reorder nodes within modules
On = zeros(n,1,'uint64');                                   %node order array
for j=1:m
    i = Om(j);
    ind = find(M==i);                                       %indices
    pos = find(Om==i);                                      %position
    
    mod_imp=[Om; sign((1:m)-pos); abs((1:m)-pos); Am(i,Om)].';
    mod_imp=sortrows(mod_imp,[3 -4]);    
    mod_imp=mod_imp(:,1).*mod_imp(:,2);
    mod_imp=[mod_imp(2:end);Om(j)];
    
	[dum ord]=sortrows(Knm(ind,:),mod_imp);                 %sort nodes by number of links to close modules
	On(ind(ord))=j*1e6+(1:Nm(i));                           %assign node order (assumes <1e6 nodes in a module)
end

[dum On]=sort(On);                                          %reorder nodes
Ar=A(On,On);                                                %reorder matrix