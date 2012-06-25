/*
    This file is part of the TOBI Interface A (TiA) library.

    Commercial Usage
    Licensees holding valid Graz University of Technology Commercial
    licenses may use this file in accordance with the Graz University
    of Technology Commercial License Agreement provided with the
    Software or, alternatively, in accordance with the terms contained in
    a written agreement between you and Graz University of Technology.

    --------------------------------------------------

    GNU Lesser General Public License Usage
    Alternatively, this file may be used under the terms of the GNU Lesser
    General Public License version 3.0 as published by the Free Software
    Foundation and appearing in the file lgpl.txt included in the
    packaging of this file.  Please review the following information to
    ensure the GNU General Public License version 3.0 requirements will be
    met: http://www.gnu.org/copyleft/lgpl.html.

    In case of GNU Lesser General Public License Usage ,the TiA library
    is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with the TiA library. If not, see <http://www.gnu.org/licenses/>.

    Copyright 2010 Graz University of Technology
    Contact: TiA@tobi-project.org
*/

/**
* @file tia_client.h
* @brief Declaration of the TiAClient class
*
**/

#ifndef TIA_CLIENT_H
#define TIA_CLIENT_H

// STL
#include <string>

#ifndef DECL_EXPORT
  #define DECL_EXPORT
#endif

namespace tia
{
// forward declarations
class DataPacket;
class TiAClientImplBase;
class SSConfig;

/**
* @class TiAClient
* @brief The TiAClient class implements a client for the SignalServer
* The client does not need much of configuration - the only thing that a user has to know is
* the server address and the control connection port of the TOBI SignalServer.
*
* Most of the functions are blocking, which means that the call will not end until
* the function has completed (time consuming) or resulted in an error.
*
* The PIMPL (Pointer to Implementation) idiom is used to achieve binary compatibility.
*
* This class is reentrant but not thread-safe.
*
* Code Example:
*
* @code
*  TiAClient client;
*
*  try
*  {
*    client.connect("127.0.0.1", 9000);
*    client.requestConfig();
*  }
*  catch (std::exception& e)
*  {
*    // error handling
*  }
*
*  SSConfig config = client.config();
*
*  // Read data packet
*  try
*  {
*    DataPacket packet;
*    // Get a packet from the server. Blocks the caller until a packet has been received.
*    client.getDataPacket(packet);
*  }
*  catch (std::exception& e)
*  {
*    // error handling
*  }
*
* @endcode
*/
class DECL_EXPORT TiAClient
{
public:
  /**
   * @brief Default Constructor
   * Call connect() to establish a connection to a TOBI SignalServer.
   * @throws
   */
  TiAClient( bool use_new_tia);
  /**
   * @brief Destructor
   */
  virtual ~TiAClient();
  /**
   * @brief Establish a connection to a TOBI SignalServer
   * The caller will be blocked until the connection has been made or an error has occurred.
   * @param address host name or the IPv4 address of the server
   * @param port the control connection port the server listens on
   * @param use_new_tia if TiA 1.0 control commands should be used (or version 0.1)
   * @throw std::ios_base::failure if the client could not
   *   connect to the server or if the client is already connected.
   */
  virtual void connect(const std::string& address, short unsigned port);

  /**
   * @brief Tells if the client is connected to the server,
   *        i.e. if the control connection has been established.
   * @return \c true if connected \c false otherwise.
   * @sa connect()
   */
  virtual bool connected() const;

  /**
   * @brief Disconnects from the server, i.e. stops receiving and closes the control connection.
   * @throw std::ios_base::failure if an occurs while disconnecting
   * \sa connect(), connected(), receiving()
   */
  virtual void disconnect();

  /**
   * @brief Requests the meta data information from the server
   * The caller will be blocked until the request has been processed or an error has occurred.
   * @throw std::ios_base::failure if the client is not connected or if an error occurred
   * sa config()
   */
  virtual void requestConfig();
  /**
   * @brief Returns the meta data information requested from the server
   * @sa requestConfig()
   */
  virtual SSConfig config() const;
  /**
   * @brief Turns the client into receiving state
   * The caller will be blocked until the client is ready to receive or an error has occurred.
   * @param use_udp_bc if \c true data will be received via UDP broadcast instead of TCP
   * @throw std::ios_base::failure if the client is not connected or if an error occurred
   * \sa stopReceiving(), receiving(), getDataPacket()
   */
  virtual void startReceiving(bool use_udp_bc);
  /**
   * @brief Tell if the client if is in receiving state
   * @return \c true if the client is in receiving state, \c false otherwise.
   * \sa startReceiving(), stopReceiving()
   */
  virtual bool receiving() const;
  /**
   * @brief Stops receiving, i.e. calling getDataPacket() will fail afterwards.
   * @throw std::ios_base::failure if the client is not connected or if an error occurred
   * \sa startReceiving(), receiving()
   */
  virtual void stopReceiving();
  /**
   * @brief Gets a packet from the server.
   * The caller will be blocked until a packet has been received or an error has occurred.
   * @param[out] packet the received data packet
   * @throw std::ios_base::failure if the client is not connected or isn't in receiving state or
   *                               if a receiving error has occurred
   * @throw std::overflow_error if an overflow occurs
   *
   */

  /**
  * @brief todo
  */
  virtual DataPacket* getEmptyDataPacket();

  virtual void getDataPacket(DataPacket& packet);
  /**
   * @brief Sets the client's data input buffer size to the given value
   * @param size the size of the input buffer in [byte]
   */
  virtual void setBufferSize(size_t size);

protected:
  TiAClientImplBase* impl_; ///< Pointer to implementation
};

} // Namespace tobiss

//-----------------------------------------------------------------------------

#endif // TIA_CLIENT_H
