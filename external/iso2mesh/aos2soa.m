function st = aos2soa(starray)
%
%    st=aos2soa(starray)
%
%    Convert an array-of-structs (AoS) to a struct-of-arrays (SoA)
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        starray: a struct array, with each subfield a simple scalar
%
%    output:
%        str: a struct, containing the same number of subfields as starray
%             with each subfield a horizontal-concatenation of the struct
%             array subfield values.
%
%    example:
%        a=struct('a',1,'b','0','c',[1,3]');
%        st=aos2soa(repmat(a,1,10))
%
%    this file is part of JSNIRF specification: https://github.com/NeuroJSON/jsnirf
%
%    License: GPLv3 or Apache 2.0, see https://github.com/NeuroJSON/jsnirf for details
%

if (nargin < 1 || ~isstruct(starray))
    error('you must give an array of struct');
end

if (numel(starray) == 1)
    st = starray;
    return
end

st = struct;
fn = fieldnames(starray);
for i = 1:length(fn)
    st.(fn{i}) = [starray(:).(fn{i})];
end
