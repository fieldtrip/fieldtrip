function as = soa2aos(starray)
%
%    as=soa2aos(starray)
%
%    Convert a struct-of-arrays (SoA) to an array-of-structs (AoS)
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        starray: a struct array, with each subfield of numeric vectors
%
%    output:
%        as: a struct array, containing the same number of subfields as starray
%             with each subfield a single scalar; the length of
%
%    example:
%        a=struct('a',[1,2],'b',[3,4]);
%        st=soa2aos(a)
%
%    this file is part of JSNIRF specification: https://github.com/NeuroJSON/jsnirf
%
%    License: GPLv3 or Apache 2.0, see https://github.com/NeuroJSON/jsnirf for details
%

if (nargin < 1 || ~isstruct(starray))
    error('you must give a struct with subfield of numeric vectors');
end

if (length(starray) > 1)
    as = starray;
    return
end

allsize = structfun(@numel, starray);
if (length(unique(allsize)) > 1)
    as = starray;
    return
end

strcell = [fieldnames(starray) struct2cell(structfun(@num2cell, starray, 'uniformoutput', 0))]';
as = struct(strcell{:});
