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
* @file data_packet_interface.h
*
* @brief Interface for the TiA datapacket
*
**/

#ifndef DATAPACKET_INTERFACE_H
#define DATAPACKET_INTERFACE_H

#include <vector>
#include <boost/cstdint.hpp>
#include <iostream>

namespace tia
{

/**
* @class DataPacket
*
* @brief Interface class for the TiA datapacket
*
*/
class DataPacket
{
  public:

    virtual ~DataPacket() {}
    virtual void reset() = 0;
    virtual void reset(void* mem) = 0;

    virtual void incPacketID() = 0;
    virtual void setPacketID(boost::uint64_t nr) = 0;
    virtual boost::uint64_t getPacketID() = 0;

     virtual void insertDataBlock(std::vector<double> v, boost::uint32_t signal_flag,
                                  boost::uint16_t blocksize, bool prepend = false)  = 0;

    virtual void setConnectionPacketNr(boost::uint64_t) = 0;
    virtual boost::uint64_t getConnectionPacketNr() = 0;

    virtual void setTimestamp() = 0;
    virtual boost::uint64_t getTimestamp() = 0;

    virtual bool hasFlag(boost::uint32_t f) = 0;

    virtual boost::uint16_t getNrOfSignalTypes() = 0;

    virtual boost::uint32_t getFlags() = 0;

    virtual std::vector<boost::uint16_t> getNrSamplesPerChannel() = 0;
    virtual std::vector<boost::uint16_t> getNrOfChannels() = 0;
    virtual std::vector<boost::uint16_t> getNrOfSamples() = 0;

    virtual boost::uint16_t getNrSamplesPerChannel(boost::uint32_t flag) = 0;
    virtual boost::uint16_t getNrOfSamples(boost::uint32_t flag) = 0;
    virtual boost::uint16_t getNrOfChannels(boost::uint32_t flag) = 0;

    virtual const std::vector<double>& getData() = 0;
    virtual std::vector<double> getSingleDataBlock(boost::uint32_t flag) = 0;

    virtual void* getRaw() = 0;
    virtual boost::uint32_t getRawMemorySize() = 0;
    virtual boost::uint32_t getRequiredRawMemorySize() = 0;
    virtual boost::uint32_t getRequiredRawMemorySize(void* mem, boost::uint32_t bytes_available) = 0;

  protected:
    DataPacket() {}

  private:
    DataPacket(const DataPacket &src) { }
    virtual DataPacket& operator=(const DataPacket &src)  {return *this;};


};

} // Namespace tia

#endif // DATAPACKET_H
