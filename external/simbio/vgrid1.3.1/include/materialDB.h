#ifndef SIMBIO_GBE_MATERIALDB_H
#define SIMBIO_GBE_MATERIALDB_H

/*
  © Copyright 2003, C&C Research Laboratories, NEC Europe Ltd.
  

    This file is part of SimBio-Vgrid.

    SimBio-Vgrid is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SimBio-Vgrid is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SimBio-Vgrid.  If not, see <http://www.gnu.org/licenses/>.


*/

/*! \file 

    Uniform access to material parameters

     $Id$

    \author Guntram Berti <berti@ccrl-nece.de>
*/


#include <string>
#include <fstream>

class materialDB_data;

/*! \brief Central data base for materials

    This class maintains a relationship between material names and id's
    used in the SimBio project, and also makes available their mechanical properties 
   (currently linear-elastic isotropic case only).

   \author Guntram Berti
*/
class materialDB {
public:
  typedef unsigned material_id_type;
  materialDB();
  materialDB(std::string   const& filenm);
  materialDB(std::ifstream      & file);
  ~materialDB();


  std::string      name(material_id_type   m) const;
  material_id_type id  (std::string const& m) const;  


  bool defined(material_id_type   m) const;
  bool defined(std::string const& m) const;  


  // isotropic elastic material constants

  //! Young's modulus  [N/mm^2]
  double E (material_id_type m) const;
  //! Poisson ratio
  double nu(material_id_type m) const;

  //! Lame constants
  double lambda(material_id_type m) const;
  double mu    (material_id_type m) const;

  //! bulk modulus (compressibility)
  double kappa (material_id_type m) const;

  //! characteristic viscoelastic relaxation time
  double t_visco (material_id_type m) const;


  //! density [g/mm^3]
  double rho   (material_id_type m) const;


  // Alternativ:
  // sauberere Trennung von Material-Modellen
  /*
  // bound to a material id
    class linear_elastic_isotropic_material {
    public:
      double  E();
      // etc.
    private:
      material_id_type  m;
      materialDB const* DB;
    };

    linear_elastic_isotropic_material material(material_id_type m) const;

   // 
   class nonlinear_elastic_isotropic_material {
   public:
     double E(stress s ???) const;
   }; 
   class linear_elastic_orthotropic_material  {
   public:
     double E(direction x ???) const;
   };
  */

private:

  //  static data d;
  materialDB_data * d;
};

#endif
