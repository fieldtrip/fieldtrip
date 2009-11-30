% function constructs matrix with all possible sets of elements; 
% [ el of SET1 , el of SET2, el of SET3 , el of SET4 ]

function [ C ] = myvariation(varargin)

% MYVARIATION produces all possible combinations of parameters from the
% input parameter lists
%
% Use as
%       [ C ] = myvariation([PARAMETER#1_list],[PARAMETER#2_list],..,[PARAMETER#n_list]) OR
%       [ C ] = myvariation({PARAMETER#1_list},{PARAMETER#2_list},..,{PARAMETER#n_list}) OR
%       
%  INPUT
%        parameter_lists as cell or numeric arrays
%
%  OUTPUT
%        array of all possible combinations of parameters from the input lists

% Pawel Herman, 2008
 
  N=1;
  n_sets = nargin;
  collective_cell = {};
  for i=1:n_sets
      
      if iscell(varargin{i})
          aux_cell = varargin{i};
          for j=1:length(aux_cell)
              N = N * length(aux_cell{j});
          end
          collective_cell = [collective_cell aux_cell];  
      elseif isnumeric((varargin{i}))
          aux_mat = varargin{i};
          N = N * length(aux_mat);
          collective_cell = [collective_cell {aux_mat}];
      else
          disp(sprintf('The parameter nr %d is neither a cell array or a numeric array',i));
          break;
      end

  end

  C=[];

  b=N;
  for k=1:length(collective_cell)  
    m = collective_cell{k};    
    b = b / length(m);        
    col=column_variation(m,b,N);
    C = [C col];  
  end
  
return;