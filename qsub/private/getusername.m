function [username] = getusername()
  if (ispc())
    username = getenv('UserName');
  else
    username = getenv('USER');
  end
  
  if (isempty(username))
    username = '<unknown>';
  end
end