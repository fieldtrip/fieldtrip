#include <winsock2.h>
#include "poll.h"

int poll (struct pollfd *p, int num, int timeout)
{
  struct timeval tv;
  fd_set read, write, except;
  int    i, n, ret;
  char buf[1024];

  FD_ZERO (&read);
  FD_ZERO (&write);
  FD_ZERO (&except);

  n = -1;
  for (i = 0; i < num; i++)
  {
    if (p[i].fd < 0)
       continue;
    if (p[i].events & POLLIN)
       FD_SET (p[i].fd, &read);
    if (p[i].events & POLLOUT)
       FD_SET (p[i].fd, &write);
    if (p[i].events & POLLERR)
       FD_SET (p[i].fd, &except);
    if (p[i].fd > n)
       n = p[i].fd;
  }

  if (n == -1)
     return (0);

  if (timeout < 0)
     ret = select (n+1, &read, &write, &except, NULL);
  else
  {
    tv.tv_sec  = timeout / 1000;
    tv.tv_usec = 1000 * (timeout % 1000);
    ret = select (n+1, &read, &write, &except, &tv);
  }

  for (i = 0; ret >= 0 && i < num; i++)
  {
    p[i].revents = 0;
    if (FD_ISSET (p[i].fd, &read))
    {

    	int j = recv(p[i].fd, buf, 1024, MSG_PEEK);
    	if(j>0)
    		p[i].revents |= POLLIN;
    	else
    		p[i].revents |= POLLHUP;
    			
    }
    if (FD_ISSET (p[i].fd, &write))
       p[i].revents |= POLLOUT;
    if (FD_ISSET (p[i].fd, &except))
       p[i].revents |= POLLERR;
  }
  return (ret); 
}
