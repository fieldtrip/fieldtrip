% FT_PREAMBLE_HELP is a helper script that display the calling-function's
% help in case the user did not specify any input argument. This can be
% used in all fieldtrip main functions that take at least a cfg input
% argument, and most also take one or multiple data structures.

if nargin==0
  stack = dbstack('-completenames');
  % stack(1) is this script
  % stack(2) is the calling ft_postamble function
  % stack(3) is the main FieldTrip function that we are interested in
  stack = stack(3);
  help(stack.name);
  % throw the error as if it happened in the original function
  msg.message     = 'This function requires one or multiple input arguments, please refer to the documentation above';
  msg.identifier  = '';
  msg.stack       = stack;
  error(msg);
end
