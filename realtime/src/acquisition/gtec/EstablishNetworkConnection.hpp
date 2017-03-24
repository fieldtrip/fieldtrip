//COPYRIGHT © 2013 G.TEC MEDICAL ENGINEERING GMBH, AUSTRIA
#ifndef ESTABLISHNETWORKCONNECTION_HPP_INCLUDED
#define ESTABLISHNETWORKCONNECTION_HPP_INCLUDED

#include <string.h>
#include <string>
#include <stdio.h>
#include <iostream>
#include <fcntl.h>

#include "gdsnetworkinterfacedemo_platform_specific.hpp"

#define NETWORKING_STARTUP_SUCCESS 0

// Initialize the networking components
//-------------------------------------------------------------------------------------
int InitNetworking();

// Free resources of the networking components if a network session ends
//-------------------------------------------------------------------------------------
void CleanupNetworking();

// Connects to the server on ip:port and sets sockefd to the associated socket file
// descriptor
//-------------------------------------------------------------------------------------
int EstablishNetworkConnection(socketd_t* socketfd, std::string ip, int port);

// Closes / disconnects the associated connection
//-------------------------------------------------------------------------------------
void CloseNetworkConnection(const socketd_t& socketfd);

// Opens a port which is listening for incoming connections
//-------------------------------------------------------------------------------------
int ListenOnNetwork(socketd_t* listenfd, int port);

// Accepts incoming connection on the socket associate with the listening socket file
// descriptor
//-------------------------------------------------------------------------------------
int AcceptOnNetwork(socketd_t listen_fd, size_t receive_buffer_size = TCP_RECV_BUF_SIZE);

#endif // ESTABLISHNETWORKCONNECTION_HPP_INCLUDED
