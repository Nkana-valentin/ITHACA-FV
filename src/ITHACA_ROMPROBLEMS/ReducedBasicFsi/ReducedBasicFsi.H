/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------
License
    This file is part of ITHACA-FV
    ITHACA-FV is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.
Class
    ReducedBasicFsi
Description
    A reduced problem for the Basic Fsi with mesh motion equations
SourceFiles
    ReducedBasicFsi.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ReducedBasicFsi class.
/// \dir
/// Directory containing the header and source files for the ReducedBasicFsi class

#ifndef ReducedBasicFsi_H
#define ReducedBasicFsi_H

#include "fvCFD.H"
#include "IOmanip.H"
//#include "ReducedProblem.H"
#include "fsiBasic.H"
#include "ITHACAutilities.H"
#include "ReducedSimpleSteadyNS.H"
#include "EigenFunctions.H"
#include "pointConstraints.H"
#include <Eigen/Eigen>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include "Modes.H"


//#include "fsiBasic.H"
#include "ITHACAPOD.H"
//#include "ReducedSimpleSteadyNS.H"
#include "ITHACAstream.H"
#include "dynamicFvMesh.H"
#include "ReducedProblem.H"
#include <chrono>
//#include<math.h>
//#include<iomanip>
//#include "pointConstraints.H"
//#include "mathematicalConstants.H"

//class fsiBasic;
class ReducedBasicFsi : public reducedSimpleSteadyNS
{
    public:
        /// Constructor
        ReducedBasicFsi();
        explicit ReducedBasicFsi(fsiBasic& problem);
        ~ReducedBasicFsi() {};

    // Variables
    /// full problem
	//fsiBasic* problem;
	//volScalarModes Pmodes;
	//volVectorModes Umodes;
	//pointVectorModes Dmodes;
    fsiBasic* problem;
    volScalarModes Pmodes;
    volVectorModes Umodes;
    pointVectorModes Dmodes;
   
    ///////// time control variables
    scalar startTime = 0.0;
    scalar finalTime = 0.0;
    scalar timeStep = 0.0;
    scalar writeEvery = timeStep;
    scalar nextWrite = 0.0;
    /// List scalar for access the centerofmass
    List<scalar> centerofmassx;
    List<scalar> centerofmassy;
    List<scalar> centerofmassz;
     /// List scalar for access the velocities of the centerofmass
    List<scalar> velx;
    List<scalar> vely;
    List<scalar> velz;
     /// List to save lift and drag forces
    List<scalar> romforcey;
    List<scalar> romforcex;
    // for saving pod interpolation  coefficients using rbf
    List<scalar> pdcoeffrbf;
    /// List of POD coefficients
    List<Eigen::MatrixXd> CoeffU;
    List<Eigen::MatrixXd> CoeffP;

    PtrList<volScalarField> PredFields;
    PtrList<volVectorField> UredFields;
    PtrList<pointVectorField> DfieldRbf;
    PtrList<pointVectorField> Dfield;
    PtrList<volScalarField> Pfield;
    PtrList<volVectorField> Ufield;



    label counter = problem->counter;

	// /// time control variables
	// scalar startTime = 0.0;
	// scalar finalTime = 0.0;
	// scalar timeStep = 0.0;
	// scalar writeEvery = timeStep;
	// scalar nextWrite = 0.0;

	// /// List scalar for access the centerOfmass
	// List<scalar> centerofmassx;
	// List<scalar> centerofmassy;
	// List<scalar> centerofmassz;

	// /// List scalar for access the velocities of the centreOfmass
	// List<scalar> velx;
	// List<scalar> vely;
	// List<scalar> velz;

	// /// List to save lift and drag forces in the reduced level
	// List<scalar> romforcey;
	// List<scalar> romforcex;
 //    /// Reduced level fields
	// PtrList<volScalarField> PredFields;
	// PtrList<volVectorField> UredFields;
	// PtrList<pointVectorField> Dfield;

	// label counter;

    // Function to perform the online phase
    void solveOnline_Pimple(scalar mu_now, int NmodesUproj, int NmodesPproj, int NmodesDproj, fileName folder = "./ITHACAoutput/Reconstruct/");
    ///
    /// @brief      It checks if the number of imposed boundary conditions is correct
    /// and set the inlet velocity equal to the given one.
    ///
    /// @param[in]  vel   Inlet velocity
    ///
    //void setOnlineVelocity(Eigen::MatrixXd vel);

    bool checkWrite(Time& timeObject);

    void PodIpointDispl(Eigen::MatrixXd coeffL2, Eigen::MatrixXd muu,  label NPdModes);

};

#endif

