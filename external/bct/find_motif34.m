function M=find_motif34(m,n)
%FIND_MOTIF34       Motif legend
%
%   Motif_matrices = find_motif34(Motif_id,Motif_class);
%   Motif_id = find_motif34(Motif_matrix);
%
%   This function returns all motif isomorphs for a given motif id and 
%   class (3 or 4). The function also returns the motif id for a given
%   motif matrix
%
%   1. Input:       Motif_id,           e.g. 1 to 13, if class is 3
%                   Motif_class,        number of nodes, 3 or 4.
%
%      Output:      Motif_matrices,     all isomorphs for the given motif
%
%   2. Input:       Motif_matrix        e.g. [0 1 0; 0 0 1; 1 0 0]
%
%      Output       Motif_id            e.g. 1 to 13, if class is 3
%
%
%Mika Rubinov, UNSW, 2007-2008

persistent M3 ID3 M4 ID4

if isscalar(m)
    if n==3
        if isempty(ID3);
            load motif34lib M3 ID3;
        end
        ind=find(ID3==m).';
        M=zeros(3,3,length(ind));
        for i=1:length(ind)
            M(:,:,i)=reshape([0 M3(ind(i),1:3) 0 ...
                M3(ind(i),4:6) 0],3,3);
        end
    elseif n==4
        if isempty(ID4);
            load motif34lib M4 ID4;
        end
        ind=find(ID4==m).';
        M=zeros(4,4,length(ind));
        for i=1:length(ind)
            M(:,:,i)=reshape([0 M4(ind(i),1:4) 0 ...
                M4(ind(i),5:8) 0 M4(ind(i),9:12) 0],4,4);
        end
    end
else
    n=size(m,1);
    M=eval(['find(motif' int2str(n) 'struct_bin(m))']);
end