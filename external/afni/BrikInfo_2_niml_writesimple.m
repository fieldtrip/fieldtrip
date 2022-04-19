function S=BrikInfo_2_niml_writesimple(Info, S)
%
%   [S]=BrikInfo_2_niml_writesimple(Info, [S])
%
%Purpose:
% Converts some fields in the BRIK header structure
% to a structure used in function afni_niml_writesimple
%
%Input Parameters:
%  Info: Brick header structure like one output by BrikInfo
%  S: If specifed add fields to this struct
%
%Output:
%  S: Struct to be used in afni_niml_writesimple
%
%

if (nargin == 1),
   S = struct();
end

if (isfield(Info,'BRICK_LABS') & ~isempty(Info.BRICK_LABS)),
   S.labels = split_string(Info.BRICK_LABS,'~');
end
if (isfield(Info,'HISTORY_NOTE') & ~isempty(Info.HISTORY_NOTE)),
   S.history = Info.HISTORY_NOTE;
end
if (isfield(Info,'BRICK_STATAUX') & ~isempty(Info.BRICK_STATAUX)),
   S.stats = Stats_Aux(Info.BRICK_STATAUX, Info.DATASET_RANK(2));
end
if (isfield(Info,'NodeIndices') & ~isempty(Info.NodeIndices)),
   S.node_indices = Info.NodeIndices;
end

return;


function s = Stats_Aux(v,nsb)
%Take Info.BRICK_STATAUX and nimlize it

%tt array is based on variable distname in niml_stat.c
%The first 'none' (there is 'none', 'none' in distname[])
%is taken out so that I can directly index into tt
%
tt = { 'none'    , 'Correl'     , 'Ttest'      , 'Ftest'    ,...
   'Zscore'   , 'Chisq'   , 'Beta'       , 'Binom'      , 'Gamma'    ,...
   'Poisson'  , 'Normal'  , 'Ftest_nonc' , 'Chisq_nonc' , 'Logistic' ,...
   'Laplace'  , 'Uniform' , 'Ttest_nonc' , 'Weibull'    , 'Chi'      ,...
   'Invgauss' , 'Extval'  , 'Pval'       , 'LogPval'    , 'Log10Pval' };

%Change stat aux format
nv = length(v);
s = repmat(cellstr('none'),1,nsb);
i = 1;
while (i<nv),
   sb = v(i); %sub-brick
   st = v(i+1); %stat type
   np = v(i+2); %number of params
   s(sb+1) = cellstr(sprintf('%s(',char(tt(st))));
   j=1;
   while(j<=np),
      if (j==np),
         s(sb+1) = cellstr(sprintf('%s%g)',char(s(sb+1)),v(i+2+j)));
      else
         s(sb+1) = cellstr(sprintf('%s%g,',char(s(sb+1)),v(i+2+j)));
      end
      j = j + 1;
   end
   i = i+3+np;
end

return;
