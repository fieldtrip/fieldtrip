/* Sample UDP server */

#include <sys/socket.h>
#include <netinet/in.h>
#include <strings.h>
#include <stdio.h>

int main(int argc, char**argv)
{
   int sockfd, n, count=0;
   struct sockaddr_in servaddr, cliaddr;
   socklen_t len;
   char mesg[2048];

   sockfd = socket(AF_INET,SOCK_DGRAM,0);

   bzero(&servaddr,sizeof(servaddr));
   servaddr.sin_family      = AF_INET;
   servaddr.sin_addr.s_addr = htonl(INADDR_ANY);
   servaddr.sin_port        = htons(55000);
   bind(sockfd,(struct sockaddr *)&servaddr,sizeof(servaddr));

   for (;;)
   {
      len = sizeof(cliaddr);
      n = recvfrom(sockfd,mesg,2048,0,(struct sockaddr *)&cliaddr,&len);
      sendto(sockfd,mesg,n,0,(struct sockaddr *)&cliaddr,sizeof(cliaddr));
      mesg[n] = 0;
      printf("Received %d bytes, count = %d\n", n, ++count);
   }
}
