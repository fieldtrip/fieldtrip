
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

    $Id: materialDB.C,v 1.10 2004/02/16 14:49:27 berti Exp $

    \author Guntram Berti <berti@ccrl-nece.de>
*/

#include "materialDB.h"

#include <vector>
#include <map>
#include <cassert>
#include <iostream>
#include <stdlib.h>

using  namespace std; 


//--------------- class materialDB_data ---------------

class materialDB_data {
  typedef materialDB::material_id_type mat_id_t;
  friend class materialDB;
public:
  materialDB_data();
  materialDB_data(std::ifstream & f);

  unsigned NumOfMaterials() const 
  { return 100; /* to be update if more materials are defined in SimBio */ }
private:
  void init_material(std::string const& name,
		     mat_id_t           m,
		     double             E_,
		     double             nu_,
		     double             kappa_,
		     double             t_visco_,
		     double             rho_
		     );


  std::map<std::string, mat_id_t> name2id;
  std::vector<std::string>        id2name;
 
  std::vector<double>  E;  // Youngs modulus [N/mm^2]
  std::vector<double>  nu; // Poisson ratio  
  
  std::vector<double>  kappa; // bulk modulus (compressbility)
  std::vector<double>  t_visco; // relaxation time
  std::vector<double>  rho;   // density     [g/mm^3]

  std::vector<bool>   defined;                    
  };



materialDB_data::materialDB_data() 
  : id2name(NumOfMaterials()+1, ""),
    E      (NumOfMaterials()+1, -1),
    nu     (NumOfMaterials()+1, -1),
    kappa  (NumOfMaterials()+1, -1),
    t_visco (NumOfMaterials()+1, -1),
    rho    (NumOfMaterials()+1, -1),
    defined(NumOfMaterials()+1, false)
{
  // initialize tables

  // IDs taken from D1.2b, page 6
  //             name    id  E          nu      kappa t_visco rho
  init_material("skull", 1,  6.50e3,    0.22,   0,    -1,     0); // from Remmler et al [1]
  init_material("scalp", 3,  1.60e-1,   0.42,   0,    -1,     0);
  init_material("brain", 5,  6.67e-2,   0.48,   0,    -1,     0);
  // to be completed
  /* 
  [1] Remmler, D. et al: Pre-surgical CT/FEA fro craniofacial distraction,
      Med Eng Phys 20 (1998), 607-619
  */
};

materialDB_data::materialDB_data(std::ifstream & f) 
  : id2name(NumOfMaterials()+1, ""),
    E      (NumOfMaterials()+1, -1),
    nu     (NumOfMaterials()+1, -1),
    kappa  (NumOfMaterials()+1, -1),
    t_visco  (NumOfMaterials()+1, -1),
    rho    (NumOfMaterials()+1, -1),
    defined(NumOfMaterials()+1, false)
{
  string    matnm; 
  mat_id_t  matid=0;
  double E_=-1, nu_=-1, kappa_=-1, t_visco_=-1, rho_ =-1;
  std:: cerr  << "materialDB: reading Materials: " << endl; 

  while(f >> matnm) {
     if(f >> matid >> E_ >> nu_ >> kappa_ >> t_visco_ >> rho_)
     {
        std:: cerr << matnm << ' ' << matid << ' ' << E_ << ' ' << nu_ << ' ' << kappa_ << ' '  << t_visco_ << ' '<< rho_ << endl;
        init_material(matnm, matid, E_, nu_, kappa_, t_visco_, rho_);
     }
     else
     {
        std:: cerr << " Something went wrong while reading the parameters of material \""
	           << matnm << "\"." << endl
		   << " The material name should be followed by an integer ID"
		   << " and 5 numbers for E, nu, kappa, t_visco and rho." << endl;
	assert(f);
     }
  }   
  std::cerr << " done materialDB: reading Materials." << endl;
}

void materialDB_data::init_material(std::string const& name,
				    mat_id_t           m,
				    double             E_,
				    double             nu_,
				    double             kappa_,
				    double             t_visco_,
				    double             rho_
				    )
{
  assert (m >= 0 && m <= NumOfMaterials());

  name2id[name] = m;
  id2name[m]    = name;

  E    [m] = E_;
  nu   [m] = nu_;
  kappa[m] = kappa_;
  t_visco[m] = t_visco_;
  rho  [m] = rho_;

  defined[m] = true;
}



//------------------ class materialDB ---------------------


materialDB::materialDB() {
  d = new materialDB_data;
}

materialDB::materialDB(std::ifstream & f) {
  d = new materialDB_data(f);
}


materialDB::materialDB(std::string const& f_nm) {
  if(f_nm == "")
    d = new materialDB_data;
  else {
    ifstream f(f_nm.c_str());
    if(! f.is_open()) {
      cerr << "ERROR: materialDB: could not open material data base file " << f_nm << endl;
      exit(2);
    }
    else {
      d = new materialDB_data(f);
      f.close();
    }
  }
}

materialDB::~materialDB() {
  delete d;
}

typedef materialDB::material_id_type mat_id;

std::string materialDB::name(mat_id m)             const { return d->id2name[m];}
mat_id      materialDB::id  (std::string const& m) const { return d->name2id[m];}

double materialDB::E (mat_id m) const { assert(defined(m));  return d->E [m];}
double materialDB::nu(mat_id m) const { assert(defined(m));  return d->nu[m];}

double materialDB::lambda(mat_id m) const { 
  double e = E(m), n = nu(m);
  return (e*n)/((1+n)*(1-2*n)); 
}
double materialDB::mu(mat_id m) const { 
  double e = E(m), n = nu(m);
  return e/(2+2*n); 
}
double materialDB::kappa(mat_id m) const {
  // For compressible material (nu < 0.5)
  // kappa = E/(3-6*nu) halds for the bulk modulus kappa.
  // But for (nearly) incompressible material the
  // bulk modulus in some element fromulations
  // for instance "mean dilatation method" plays
  // the role of a "penalty term" and must be set
  // differently (for instance kappa = 10^4 * mu
  // for the "mean dilatation method".
  // Therefore kappa is read from the material database,
  // but the user has to make sure, that its value is consistent
  // with its use !!!!! 
  
  assert(defined(m)); return d->kappa[m]; 
}

double materialDB::t_visco  (mat_id m) const { assert(defined(m));  return d->t_visco  [m];}

double materialDB::rho  (mat_id m) const { assert(defined(m));  return d->rho  [m];}

bool materialDB::defined(mat_id m) const 
{ return ((m >= 0) && (m <= d->NumOfMaterials()) && d->defined[m]);}

