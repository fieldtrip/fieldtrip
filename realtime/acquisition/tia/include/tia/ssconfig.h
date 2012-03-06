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
* @file ssconfig.h
* @brief This file includes a class storing the TiA config for the subject and the transmitted signals.
**/

#ifndef SSCONFIG_BASE_H
#define SSCONFIG_BASE_H

// local
#include "tia/ss_meta_info.h"

namespace tia
{
//-----------------------------------------------------------------------------
/**
* @class SSConfig
*
* @brief This class holds information regarding the subject and the transmitted signals.
*
* @todo Rename class to TiAConfig
*/
class SSConfig
{
  public:
  /**
  * @brief Element holding subject specific information.
  */
  SubjectInfo subject_info;
  /**
  * @brief Element holding information of transmitted signals.
  */
  SignalInfo  signal_info;
};

} // Namespace tobiss

//-----------------------------------------------------------------------------

#endif // SSCONFIG_BASE_H
