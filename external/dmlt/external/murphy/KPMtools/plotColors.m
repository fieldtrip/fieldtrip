
colors =  ['r' 'b' 'k' 'g' 'c' 'y' 'm' ...
	   'r' 'b' 'k' 'g' 'c' 'y' 'm'];
symbols = ['o' 'x' '+' '>' '<' '^' 'v' ...
	   '*' 'p' 'h' 's' 'd' 'o' 'x'];
for i=1:length(colors)
  styles{i} = sprintf('-%s%s', colors(i), symbols(i));
end

