  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [seq] = fun_shift(seq, step, dim)
  % step < 0    -----  Left or up 
  % step > 0    -----  Right or Down
  % Dim =1      -----  Row
  % Dim =2      -----  Col
  [row col]=size(seq);
  
  step = -1*round(step);
  
  if (dim==2) && (abs(step)>=col)
      seq=zeros(row, col);
      return;
  end;
  
  if (dim~=2) && (abs(step) >= row)
        seq=zeros(row, col);
     return;
  end;
  
  if step>0
      if dim ==2
          seq= [seq(:, step+1:end) zeros(row, step)];
      else 
          seq= [seq(step+1:end,:); zeros(step, col)];
      end;
  end;
  
  if step <0
      if dim ==2
          seq= [zeros(row, abs(step)) seq(:, 1:end+step) ];
      else 
          seq= [zeros(abs(step), col); seq(1:end+step,:) ];
      end;
  end;