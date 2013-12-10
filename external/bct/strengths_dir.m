function [is,os,str] = strengths_dir(CIJ)
%STRENGTHS_DIR      Instrength and outstrength
%
%   [is,os,str] = strengths_dir(CIJ);
%
%   Node strength is the sum of weights of links connected to the node. The
%   instrength is the sum of inward link weights and the outstrength is the
%   sum of outward link weights.
%
%   Input:      CIJ,    directed weighted connection matrix
%
%   Output:     is,     node instrength
%               os,     node outstrength
%               str,    node strength (instrength + outstrength)
%
%   Notes:  Inputs are assumed to be on the columns of the CIJ matrix.
%
%
%   Olaf Sporns, Indiana University, 2002/2006/2008


% compute strengths
is = sum(CIJ,1);    % instrength = column sum of CIJ
os = sum(CIJ,2)';   % outstrength = row sum of CIJ
str = is+os;        % strength = instrength+outstrength


