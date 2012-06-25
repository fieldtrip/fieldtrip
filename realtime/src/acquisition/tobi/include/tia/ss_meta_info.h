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
* @file ss_meta_info.h
* @brief This file includes classes used to store meta information used in TiA.
* @todo Rename file to TiA_meta_information
**/

#ifndef SERVER_META_INFO_H
#define SERVER_META_INFO_H

#include <boost/cstdint.hpp>

// STL
#include <string>
#include <map>
#include <vector>

namespace tia
{

//-----------------------------------------------------------------------------

/**
 * @class SubjectInfo
 *
 * @brief This class is used to store subject specific information.
 *
 * Subjects with obvious functionality are not documented.
 * @todo Adjust class to use enums properly (e.g. sex).
 * @todo Check if mandatory settings are done (now done in XML parser).
 */

class SubjectInfo
{
  public:
    enum Sex {
      Male   = 1,
      Female = 2
    };

    enum Handedness {
      RightHanded = 1,
      LeftHanded = 2
    };

    /**
    * @brief Enum for binary information.
    */
    enum ShortInfoType{
      Glasses = 1,
      Smoking
    };

    enum ShortInfoValue {
      Unknown = 0,
      No,
      Yes
    };

    typedef std::map<ShortInfoType, ShortInfoValue> ShortInfoMap;

  public:
    /**
    * @brief Get the subjects identification code.
    */
    std::string id() const { return id_; }
    /**
    * @brief Set the subjects identification code.
    */
    void setId(const std::string& id) { id_ = id; }

    std::string firstName() const { return firstname_; }
    void setFirstName(const std::string& name) { firstname_ = name; }

    std::string surname() const { return surname_; }
    void setSurname(const std::string& name) { surname_ = name; }

    std::string birthday() const { return birthday_; }
    void setBirthday(const std::string& birthday) { birthday_ = birthday; }

    int sex() const { return sex_; }
    void setSex(int sex) { sex_ = sex; }

    int handedness() const { return handedness_; }
    void setHandedness(int handedness) { handedness_ = handedness; }

    void setMedication(std::string& medication) { medication_ =  medication; }
    std::string medication() const { return medication_; }

    /**
     * @brief Set a single short info.
     */
    void setShortInfo(ShortInfoType info, ShortInfoValue value)
    {
      short_infos_[info] = value;
    }

    /**
     * @brief Get a specific short info value.
     */
    ShortInfoValue shortInfo(ShortInfoType info) const
    {
      ShortInfoMap::const_iterator it = short_infos_.find(info);
      if (it == short_infos_.end()) return Unknown;
      return (*it).second;
    }

    ShortInfoMap& shortInfoMap() { return short_infos_; }

    const ShortInfoMap& shortInfoMap() const { return short_infos_; }

  private:
    std::string  id_;           ///< The subjects identification code.
    std::string  firstname_;    ///< Voluntary.
    std::string  surname_;      ///< Voluntary.
    std::string  birthday_;     ///< Mandatory
    std::string  medication_;   ///< Mandatory
    int handedness_;            ///< Mandatory
    int sex_;                   ///< Mandatory
    ShortInfoMap short_infos_;  ///<
};

//-----------------------------------------------------------------------------

/**
 * @class Channel
 *
 * @brief Information for a channel.
 *
 * This class is available to store channel specific information. Up to now
 * only the channel name is supported. It will be extended in near future with
 * physical and digital range, impedance and related hardware ID.
 *
 * @todo Write (or rename to) more appropriate method-names(e.g. id_--> name_)
 */
class Channel
{
  public:
    /**
     * @brief Get the channel name (e.g. C3).
     */
    std::string id() const { return id_; }
    /**
     * @brief Set the channel name.
     */
    void setId(const std::string& id) { id_ = id; }

  private:
    std::string id_;    ///<
};

//-----------------------------------------------------------------------------

/**
 * @class Signal
 *
 * @brief This class is written to store information for a specific signal type.
 *
 * This class is designed to hold information which is consistent
 * within a single signal type like sampling rate or block size.
 *
 * @todo Find a way to be consistent with signal types and signal type names within TiA.
 */
class Signal
{
  public:
      /**
       * @brief Set the name of the signal type.
       */
      void setType(const std::string& type) { type_ = type; }

      /**
       * @brief Get the name of the stored signal type.
       */
      std::string type() const { return type_; }

      /**
       * @brief Get the blocksize of the stored signal type.
       */
      boost::uint16_t blockSize() const { return block_size_; }

      /**
       * @brief Set the blocksize of the stored signal type.
       */
      void setBlockSize(boost::uint16_t block_size) { block_size_ = block_size; }

      /**
       * @brief Get the sampling rate of the stored signal type.
       */
      boost::uint16_t samplingRate() const { return sampling_rate_; }

      /**
       * @brief Set the sampling rate of the stored signal type.
       */
      void setSamplingRate(boost::uint16_t sampling_rate) { sampling_rate_ = sampling_rate; }

      /**
       * @brief Get a vector holding specific information for every channel of this signal type.
       */
      const std::vector<Channel>& channels() const { return channels_; }

      std::vector<Channel>& channels() { return channels_; }

  private:
    std::string           type_;
    boost::uint16_t       block_size_;
    boost::uint16_t       sampling_rate_;
    std::vector<Channel>  channels_;
};

//-----------------------------------------------------------------------------

/**
 * @class SignalInfo
 *
 * @brief This class stores information for all transmitted signals.
 *
 * This class stores information for all signal transmitted via TiA.
 * It also stores information concerning the master device, as this
 * device is responsible for the rate, data packets are delivered with.
 */
class SignalInfo
{
  public:
    typedef std::map<std::string, Signal> SignalMap;

  public:
    /**
     * @brief Get the sampling rate of the device defined as master.
     */
    boost::uint16_t masterSamplingRate() const { return master_sampling_rate_; }

    /**
     * @brief Set the sampling rate of the device defined as master.
     */
    void setMasterSamplingRate(boost::uint16_t rate) { master_sampling_rate_ = rate; }

    boost::uint16_t masterBlockSize() const { return master_block_size_; }

    void setMasterBlockSize(boost::uint16_t size) { master_block_size_ = size; }

    const SignalMap& signals() const { return signals_; }

    SignalMap& signals() { return signals_; }

  private:
    boost::uint16_t  master_block_size_;      ///<
    boost::uint16_t  master_sampling_rate_;   ///<
    SignalMap signals_;                       ///<
};

} // Namespace tobiss

//-----------------------------------------------------------------------------

#endif // UDPSERVER_H
