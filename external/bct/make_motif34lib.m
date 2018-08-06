function make_motif34lib
%MAKE_MOTIF34LIB        Auxiliary motif library function
%
%   make_motif34lib;
%
%   This function generates the motif34lib.mat library required for all
%   other motif computations.
%
%
%   Mika Rubinov, UNSW, 2007-2010

%#ok<*ASGLU>

[M3,M3n,ID3,N3]=motif3generate; 
[M4,M4n,ID4,N4]=motif4generate;
save motif34lib;

function [M,Mn,ID,N]=motif3generate
n=0;
M=false(54,6);                  %isomorphs
CL=zeros(54,6,'uint8');         %canonical labels (predecessors of IDs)
cl=zeros(1,6,'uint8');
for i=0:2^6-1                   %loop through all subgraphs
    m=dec2bin(i);
    m=[num2str(zeros(1,6-length(m)), '%d') m];  %#ok<AGROW>
    G=str2num ([ ...
        '0'   ' '  m(3)  ' '  m(5) ;
        m(1)  ' '  '0'   ' '  m(6) ;
        m(2)  ' '  m(4)  ' '  '0'   ]);         %#ok<ST2NM>
    Ko=sum(G,2);
    Ki=sum(G,1).';
    if all(Ko|Ki),              %if subgraph weakly-connected
        n=n+1;
        cl(:)=sortrows([Ko Ki]).';
        CL(n,:)=cl;             %assign motif label to isomorph
        M(n,:)=G([2:4 6:8]);
    end
end
[u1,u2,ID]=unique(CL,'rows');   %convert CLs into motif IDs

%convert IDs into Sporns & Kotter classification
id_mika=  [1  3  4  6  7  8  11];
id_olaf= -[3  6  1 11  4  7   8];
for id=1:length(id_mika)
    ID(ID==id_mika(id))=id_olaf(id);
end
ID=abs(ID);

[X,ind]=sortrows(ID);
ID=ID(ind,:);               %sort IDs
M=M(ind,:);                 %sort isomorphs
N=sum(M,2);                 %number of edges
Mn=uint32(sum(repmat(10.^(5:-1:0),size(M,1),1).*M,2));  %M as a single number

function [M,Mn,ID,N]=motif4generate
n=0;
M=false(3834,12);               %isomorphs
CL=zeros(3834,16,'uint8');      %canonical labels (predecessors of IDs)
cl=zeros(1,16,'uint8');
for i=0:2^12-1                  %loop through all subgraphs
    m=dec2bin(i);
    m=[num2str(zeros(1,12-length(m)), '%d') m];     %#ok<AGROW>
    G=str2num ([ ...
        '0'   ' '  m(4)  ' '  m(7)  ' '  m(10) ;
        m(1)  ' '  '0'   ' '  m(8)  ' '  m(11) ;
        m(2)  ' '  m(5)  ' '  '0'   ' '  m(12) ;
        m(3)  ' '  m(6)  ' '  m(9)  ' '  '0'    ]); %#ok<ST2NM>
    Gs=G+G.';
    v=Gs(1,:);
    for j=1:2,
        v=any(Gs(v~=0,:),1)+v;
    end
    if v                        %if subgraph weakly connected
        n=n+1;
        G2=(G*G)~=0;
        Ko=sum(G,2);
        Ki=sum(G,1).';
        Ko2=sum(G2,2);
        Ki2=sum(G2,1).';
        cl(:)=sortrows([Ki Ko Ki2 Ko2]).';
        CL(n,:)=cl;             %assign motif label to isomorph
        M(n,:)=G([2:5 7:10 12:15]);
    end
end
[u1,u2,ID]=unique(CL,'rows');   %convert CLs into motif IDs
[X,ind]=sortrows(ID);   
ID=ID(ind,:);                   %sort IDs
M=M(ind,:);                     %sort isomorphs
N=sum(M,2);                     %number of edges
Mn=uint64(sum(repmat(10.^(11:-1:0),size(M,1),1).*M,2)); %M as a single number