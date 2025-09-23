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
Description
Tutorial of compressible and unsteady flow around a moving airfoil
\*---------------------------------------------------------------------------*/

#include "CompressibleUnSteadyRhoPimple.H"
#include "ITHACAPOD.H"
#include "ITHACAstream.H"
#include "Foam2Eigen.H"
#include "DEIM.H"
//#include "ReducedProblem.H"
#include <chrono>


class tutorial26: public CompressibleUnSteadyRhoPimple
{
public:
        tutorial26(int argc, char* argv[])
        : CompressibleUnSteadyRhoPimple(argc, argv), 
        U(_U()), p(_p()), E(_E())
        //,pd(_pointDisplacement())
    {
        //point0 = meshPtr().points();
    }
    /// Velocity field
    volVectorField& U;
    /// Pressure field
    volScalarField& p;
    /// Energy field
    volScalarField& E;
    /// Hyper-reduced objects
    //autoPtr<MyDEIM<volScalarField>> DeimP, DeimE;
    //autoPtr<MyDEIM<volVectorField>> DeimU;
    ///grid nodes field
    //pointVectorField& pd;
    /// Initial coordinates of the grid points
    void offlineSolve(word folder="./ITHACAoutput/Offline/")
    {
        //List<scalar> mu_now(1);

        if (offline)
        {
            ITHACAstream::read_fields(Ufield, U, folder);
            ITHACAstream::read_fields(Pfield, p, folder);
            ITHACAstream::read_fields(Efield, E, folder);
        }
        else
        {
            truthSolve(folder);
            //restart();

        }

    }

}; 
    
/*----------------------------------------------------------------------------------------------------------*\
                               Starting the MAIN
\*-----------------------------------------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial26 example(argc, argv);
    tutorial26 online(argc, argv);
    
    std::clock_t startOff;
    double durationOff;
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance
                             (
                                example.meshPtr(),
                                example._runTime()
                              );
    //Info <<  para->runTime.time().constant() << endl;                         
    
    int NmodesUout  =  readInt(para->ITHACAdict->lookup("NmodesUout"));
    int NmodesPout  =  readInt(para->ITHACAdict->lookup("NmodesPout"));
    int NmodesEout  =  readInt(para->ITHACAdict->lookup("NmodesEout"));

    int NmodesUproj  = readInt(para->ITHACAdict->lookup("NmodesUproj"));
    int NmodesPproj  = readInt(para->ITHACAdict->lookup("NmodesPproj"));
    int NmodesEproj  = readInt(para->ITHACAdict->lookup("NmodesEproj"));
    // word filename("./par");
    // example.mu = ITHACAstream::readMatrix(filename);
    // Time parameters: We can use Ioodictionnary to access time parameters
    example.startTime  = 0;
    example.finalTime  = 0.15;
    example.timeStep   = 2e-06; 
    example.writeEvery = 4e-04;

    // //Perform the offline solve
    startOff= std::clock();
    example.offlineSolve();
    //exit(0);
    durationOff = (std::clock() - startOff);
    std::cout << "The Offline phase  duration  is  equal  to " << durationOff << std::endl;
 /*
    if(example.podex==0 )
    {
       ITHACAPOD::getModes(example.Ufield, online.Umodes, example._U().name(),
                    example.podex, 0, 0, NmodesUout);
       ITHACAPOD::getModes(example.Pfield, online.Pmodes, example._p().name(),
                            example.podex, 0, 0,NmodesPout);
       ITHACAPOD::getModes(example.Efield, online.Emodes, example.E().name(),
                            example.podex, 0, 0, NmodesEout);
    }
    else
    {
      ITHACAstream::read_fields(online.Umodes, example._U(), "./ITHACAoutput/POD/");
      ITHACAstream::read_fields(online.Pmodes, example._p(), "./ITHACAoutput/POD/");
      ITHACAstream::read_fields(online.Emodes, example._E(), "./ITHACAoutput/POD/");

    }*/
    
    exit(0);
}


