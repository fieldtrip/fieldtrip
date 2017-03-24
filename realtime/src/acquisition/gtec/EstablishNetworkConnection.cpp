//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#include "EstablishNetworkConnection.hpp"

// these are from FieldTrip/realtime/buffer
#include "platform.h"
#include "compiler.h"
#include "platform_includes.h"

#if !defined(PLATFORM_WINDOWS)

int InitNetworking()
{
    return 0;
}

void CleanupNetworking()
{
}

#elif defined(PLATFORM_WINDOWS)

int InitNetworking()
{
    WSADATA impinfo;
    int ret = WSAStartup( MAKEWORD(2,2), &impinfo );
    if ( ret != NETWORKING_STARTUP_SUCCESS )
    std::cerr << "ERROR: could not start winsock, error = " << ret << std::endl;

    return ret;
}

void CleanupNetworking()
{
    // clean up winsock lib
    WSACleanup();
}

#endif

int EstablishNetworkConnection(socketd_t* socketfd, std::string ip, int port)
{
    struct sockaddr_in gds_address;
    *socketfd = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);

    if (*socketfd < 0)
	{
        std::cerr << "ERROR: could not open socket " << port << std::endl;
        return (int) *socketfd;
    }

	//disable TCP's delayed acknowledgement feature (i.e. Nagle's algorithm; this is required to avoid a mess up of the general timing during data acquisition)
	int tcp_nodelay = 1;
	int ret = setsockopt(*socketfd, IPPROTO_TCP, TCP_NODELAY, (char*) &tcp_nodelay, sizeof(int));

	if (ret == SOCKET_ERROR)
        std::cerr << "ERROR: Could not disable TCP's delayed acknowledgement feature" << std::endl;

    gds_address.sin_family = AF_INET;
    gds_address.sin_addr.s_addr = inet_addr(ip.c_str());
    gds_address.sin_port = htons(port);

    ret = connect(*socketfd, reinterpret_cast<const struct sockaddr*>(&gds_address), sizeof(gds_address));
    if (ret == SOCKET_ERROR)
        std::cerr << "ERROR: could not connect to " << ip << ":" << port << std::endl;

    return ret;
}

void CloseNetworkConnection(const socketd_t &socketfd)
{
    MKR_CLOSE_SOCKET(socketfd);
}

int ListenOnNetwork(socketd_t* listenfd, int port)
{
    struct sockaddr_in local_address;
    *listenfd = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
    if (*listenfd == SOCKET_ERROR)
	{
        std::cerr << "ERROR: could not open a socket for listening" << std::endl;
        return (int) *listenfd;
    }

	int ret;

	//disable TCP's delayed acknowledgement feature (i.e. Nagle's algorithm; this is required to avoid a mess up of the general timing during data acquisition)
	int tcp_nodelay = 1;
	ret = setsockopt(*listenfd, IPPROTO_TCP, TCP_NODELAY, (char*) &tcp_nodelay, sizeof(int));

	if (ret == SOCKET_ERROR)
        std::cerr << "ERROR: Could not disable TCP's delayed acknowledgement feature" << std::endl;

#if defined linux
    /* Enable address reuse */
    int on = 1;
    ret = setsockopt(*listenfd, SOL_SOCKET, SO_REUSEADDR, &on, sizeof(on));
    if (ret < 0)
	{
        std::cerr << "ERROR: could not set socket option: REUSEADDR" << std::endl;
        return ret;
    }
#endif //DEFINED_LINUX
    memset((char*) &local_address, 0, sizeof(local_address));
    local_address.sin_family = AF_INET;
    local_address.sin_addr.s_addr = INADDR_ANY;
    local_address.sin_port = htons(port);

    ret = bind(*listenfd, (struct sockaddr *) &local_address, sizeof(local_address));
    if (ret < 0)
	{
        std::cerr << "ERROR: could not bind listening socket" << std::endl;
        return ret;
    }

    listen(*listenfd, 1);
    return (int) *listenfd;
}

int AcceptOnNetwork(socketd_t listen_fd, size_t receive_buffer_size)
{
    struct sockaddr_in client_addr;
    socklen_t client_addr_length = 0;
	int tcp_nodelay = 1;

    socketd_t socketfd = MKR_ACCEPT(listen_fd, (struct sockaddr *) &client_addr, &client_addr_length);
    if (socketfd == SOCKET_ERROR)
        std::cerr << "ERROR: on accept" << std::endl;

    int ret = setsockopt(socketfd, SOL_SOCKET, SO_RCVBUF, (char*) &receive_buffer_size, sizeof(receive_buffer_size));
    if (ret == SOCKET_ERROR)
        std::cerr << "ERROR: Could not set the receive buffer size" << std::endl;

	ret = setsockopt(socketfd, IPPROTO_TCP, TCP_NODELAY, (char*) &tcp_nodelay, sizeof(int));
	if (ret == SOCKET_ERROR)
        std::cerr << "ERROR: Could not disable TCP's delayed acknowledgement feature" << std::endl;

    return (int) socketfd;
}
