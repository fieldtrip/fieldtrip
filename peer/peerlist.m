function list = peerlist

% PEERLIST returns the list of peers and/or displ

list = peer('peerlist');

if nargout==0
  % display the hosts on screen, sort them by the hostid
  [dum, indx] = sort([list.hostid]);
  list = list(indx);
  for i=1:numel(list)
    switch list(i).hoststatus
      case 0
        status = 'zombie';
      case 1
        status = 'slave ';
      case 2
        status = 'master';
      otherwise
        error('unknown status');
    end
    fprintf('%s at %s@%s:%d, group = %s, hostid = %d\n', status, list(i).user, list(i).hostname, list(i).hostport, list(i).group, list(i).hostid);
  end
  clear list
end
