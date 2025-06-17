//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
//      Official webSite: https://code-mphi.github.io/ECOGEN/
//
//  This file is part of ECOGEN.
//
//  ECOGEN is the legal property of its developers, whose names 
//  are listed in the copyright file included with this source 
//  distribution.
//
//  ECOGEN is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published 
//  by the Free Software Foundation, either version 3 of the License, 
//  or (at your option) any later version.
//  
//  ECOGEN is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with ECOGEN (file LICENSE).  
//  If not, see <http://www.gnu.org/licenses/>.

#ifndef TOOLS_H
#define TOOLS_H

#include "Eos/Eos.h"

//#define DEBUG

//! \brief     Enumeration for the primitive-variable type (usefull for second order, slopes, etc.)
enum Prim { vecPhases, vecPhasesO2, vecSlopes, restart };

//! \brief     Enumeration for the axes (X, Y and Z for the axes in the x-, y- and z-direction)
enum Axis { X, Y, Z };

//! \brief     Enumeration for the type of output (GNU, XML)
enum TypeOutput { GNU, XML };

//! \brief     Enumeration for the type of mesh (REC: rectilinear, UNS: unstructured, AMR: adaptative mesh refinement)
enum TypeM { REC, UNS, AMR };

//! \brief     Enumeration for the type of data (FLOAT, DOUBLE, INT, CHAR)
enum TypeData { FLOAT, DOUBLE, INT, CHAR };

//! \brief     Enumeration for the type of geometric object (VERTEX, LINE, PLAN)
enum TypeGO { VERTEX, LINE, PLAN };

//! \brief     Enumeration for the type of boundary
enum TypeBC { NONREFLECTING = 1, WALL = 2, OUTLETPRESSURE = 3, INLETINJSTAGSTATE = 4, INLETTANK = 5, SYMMETRY = 6, INLETINJTEMP = 7, NULLFLUX = 8, OUTLETMASSFLOW = 9, PISTON = 10 };

//! \brief     Enumeration for the variable to extract on a boundary
enum VarBoundary { p, rho, velU, velV, velW, SIZE = (velW+1) };

//! \brief     Enumeration for the heat type of wall boundary (ADIABATIC, IMPOSEDTEMP, IMPOSEDFLUX)
enum TypeBCHeat { ADIABATIC = 0, IMPOSEDTEMP = 1, IMPOSEDFLUX = 2};

//! \brief     Enumeration for the type of relaxation (P, PT, PTMu)
enum TypeRelax { P = 1, PT = 2, PTMU = 3, GPR = 4 };

//! \brief     Enumeration for the gradient method (Finite-Difference (FD), Green-Gauss (GG))
enum TypeGrad { FD, GG };

//! \brief     Enumeration for the phase index of a liquid/vapor couple
enum PhaseType { LIQ = 0, VAP = 1 };

//! \brief     Enumeration for the type of second order limiter (NS only)
enum LimiterType { NONE = 0, MINMOD = 1, SUPERBEE = 2 };

//! \brief     Enumeration for the flow variables
enum class Variable { transport, pressure, density, alpha, velocityMag, velocityU, velocityV, velocityW, temperature, QPA, lambda,
                      cobaseXX, cobaseXY, cobaseXZ, cobaseYX, cobaseYY, cobaseYZ, cobaseZX, cobaseZY, cobaseZZ };

//! \brief     Template for the type of the mesh container (std::list for now, but may change to something else if wanted)
template<class Type>
using TypeMeshContainer=std::vector<Type>;

//! \brief     Template for type-safe sign function
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

//! \class     Tools
//! \brief     Class for tools
class Tools
{
  public:
    //! \brief     Generic model constructor
    //! \param     numbPhases         number of phases
    //! \param     numbSolids         number of solid phases
    //! \param     numberTransports   number of additional transport equations
    Tools(const int& numbPhases, const int& numbSolids, const int& numbTransports);
    ~Tools();

    //! \brief     Modify the string of characters to uppercase it
    //! \param     string               string of characters
    static void uppercase(std::string& string);
    //! \brief     Modify the string of characters to lowercase it
    //! \param     string               string of characters
	  static void lowercase(std::string& string);

    //! \brief     Swap two numbers
    //! \param     double   1st number to swap
    //! \param     double   2nd number to swap
    static void swap(double &a, double &b);

    //! \brief     Return a non-zero value of a float
    //! \param     double   initial float
    double returnNonZeroValue(double a);

    double* ak;
    double* Yk;
    double* rhok;
    double* pk;
    double* ek;
    double* Ek;
    double* akS;
    double* rhokS;
    double* rhokStar;
    double* pkStar;
    double* ekStar;
    double* EkStar;
    double* vkStar;
    double* YkStar;
    double* Deltapk;           //!< Pressure differences, one for each phase
    double* zk;                //!< Acoustic impedance of each phase
    double* rho_cIksquare;     //!< Density times interface sound speed square, one for each phase
    Eos** eos;
    double* Hk0;
    double* Yk0;
    double* compactionPk;      //!< Compaction pressure, one for each phase
    double* dlambda;           //!< Plastic compaction term, one for each phase
    Tensor* dplast;            //!< Plastic shear term, one for each phase
    bool* alphaNull;           //!< Parameter to know if we consider the volume fraction as null or not
    bool* relaxSolidPlast;     //!< Parameter to know if we need the solid plastic relaxation

    static double uselessDouble;
    double physicalTime;       //!< Current physical time
};

extern Tools *TB;
extern int numberPhases;
extern int numberSolids;
extern int numberTransports;


#endif // TOOLS_H
