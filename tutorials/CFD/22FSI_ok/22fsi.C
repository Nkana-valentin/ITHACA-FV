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
    Example of an unsteady NS Reduction Problem
SourceFiles
    22fsi.C
\*---------------------------------------------------------------------------*/

#include "fsiBasic.H"
#include "ITHACAPOD.H"
#include "ReducedSimpleSteadyNS.H"
#include "ITHACAstream.H"
#include "dynamicFvMesh.H"
#include "ReducedProblem.H"
#include <chrono>
#include<math.h>
#include<iomanip>

class tutorial22: public fsiBasic
{
public:
    explicit tutorial22(int argc, char* argv[])
        : fsiBasic(argc, argv), U(_U()), p(_p()), phi(_phi())
    {
        //point0 = meshPtr().points();
    }

    // Fields To Perform
    /// Velocity field
    volVectorField& U;
    /// Pressure field
    volScalarField& p;
    /// flux field
    surfaceScalarField& phi;
    /// Initial coordinates of the grid points
        //vectorField point0;
    void offlineSolve()
    {
        Vector<double> inl(1, 0, 0);
        List<scalar> mu_now(1);

        // if (offline)
        // {
        //     ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
        //     ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
        //      // mu_samples =
        //      //    ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
        // }
        // else
        //{
        //for (label i = 0; i < mu.cols(); i++)
        //{
        //inl[0] = mu(0, i);
        mu_now[0] = mu(0, 0); //mu.cols()=50
        //mu_now[0] = mu(0, i);
        //std::cout << "////////////////////////////mu is :///////////////////"<< mu.cols() <<  std::endl;
        //std::cout << "////////////////////////////mu_now[0] is :///////////////////"<< mu_now[0] <<  std::endl;

        //assignBC(U, BCind, inl);
        //assignIF(U, inl);
        //change_viscosity(mu(0, i));
        //meshPtr().movePoints(point0);

        truthSolve3(mu_now);
        //restart();

        //}
        //}
    }
};


class reducedBasicFsi: public reducedSimpleSteadyNS
{
public:
    explicit reducedBasicFsi(tutorial22& FoamPb): problem(&FoamPb)
        
    {
            for (int i = 0; i < problem->Umodes.size(); i++)
		    {
		        Umodes.append((problem->Umodes.toPtrList()[i]).clone());
		    }

		    for (int i = 0; i < problem->Umodes.size(); i++)
		    {
		        Pmodes.append((problem->Pmodes.toPtrList()[i]).clone());
		    }

        std::cout << "################ ctor of reducedBasicFsi ##################" << std::endl;
    }

    tutorial22* problem;
    volScalarModes Pmodes;
    volVectorModes Umodes;
   
    ///////// time control variables
    scalar startTime;
    scalar finalTime;
    scalar timeStep;
    scalar writeEvery = timeStep;
    scalar nextWrite;

   
    /// List scalar for access the centerofmass
    List<scalar> centerofmassx;
    List<scalar> centerofmassy;
    List<scalar> centerofmassz;

     /// List scalar for access the velocities of the centerofmass
    List<scalar> velx;
    List<scalar> vely;
    List<scalar> velz;

   PtrList<volScalarField> PredFields;
   PtrList<volVectorField> UredFields;
   PtrList<pointVectorField> DFields;

    label counter = 1;

    void solveOnline_Pimple(scalar mu_now, int NmodesUproj, int NmodesPproj,fileName folder = "./ITHACAoutput/Reconstruct/")
    {
        
        problem->restart();
        // Info << problem->_U()<< endl;
        // Info << problem->_p() << endl;
        Eigen::VectorXd uresidual = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::VectorXd presidual = Eigen::VectorXd::Zero(NmodesPproj);

        Eigen::MatrixXd a = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::MatrixXd a0 = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::MatrixXd b = Eigen::VectorXd::Zero(NmodesPproj);
        Eigen::MatrixXd bold =b;

       
        Time& runTime = problem->_runTime();
        dynamicFvMesh& mesh = problem->meshPtr();
//#include "initContinuityErrs.H"
        fv::options& fvOptions = problem->_fvOptions();
        volScalarField& p = problem->_p();
        volVectorField& U = problem->_U();
        surfaceScalarField& phi = problem->_phi();
        pimpleControl& pimple = problem->_pimple();
        IOMRFZoneList& MRF = problem->_MRF();
        singlePhaseTransportModel& laminarTransport = problem->_laminarTransport();
        autoPtr<incompressible::turbulenceModel> turbulence = problem->turbulence;

        // Declare modal coefficients for velocity and pressure
        //a = ITHACAutilities::getCoeffs(U, Umodes, NmodesUproj, true);
        //a(0) = a0(0);
        //b = ITHACAutilities::getCoeffs(p, Pmodes, NmodesPproj, true);

        instantList Times = runTime.times();
        runTime.setEndTime(finalTime);
        runTime.setTime(Times[1], 1);
        runTime.setDeltaT(timeStep);
        nextWrite = startTime;
        label pRefCell = 0;
        scalar pRefValue = 0.0;
        scalar  cumulativeContErr =0.0;
        //P.rename("p");

        IOdictionary   dynamicMeshDict
	    (
	         IOobject
	         (
	            "dynamicMeshDict",
	            mesh.time().constant(),
	            mesh,
	            IOobject::MUST_READ,
	            IOobject::NO_WRITE,
	            false
	         )
	    );
    /// construct a sixDoFRigidBodyMotionSolver object
    sixDoFRigidBodyMotionSolver sDRBMS(mesh, dynamicMeshDict);
    std::cerr << "/////////////////////"<< sDRBMS.motion().mass() << "///////////////////" << std::endl;

   
#include "addCheckCaseOptions.H"
#include "createDyMControls.H"
#include "createUfIfPresent.H"
#include "CourantNo.H"

        //turbulence->validate();
        //Umodes.reconstruct(U, a, "U");
        //Pmodes.reconstruct(p, b, "p");

        // ITHACAstream::exportSolution(U, name(counter), folder);
        // ITHACAstream::exportSolution(p, name(counter), folder);
        // //ITHACAstream::exportSolution(sDRBMS.pointDisplacement(), name(counter), folder);
        // ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
        // std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
        // UredFields.append(U.clone());
        // PredFields.append(p.clone());
        // //DFields.append(sDRBMS.pointDisplacement().clone());
        // counter++;
        // nextWrite += writeEvery;

        // PIMPLE algorithm starts here
        Info<< "\nStarting time loop\n" << endl;
        while (runTime.run())
        {
            std::cerr << "File: 22fsi.C, Line: 230"<< std::endl;
#include "readDyMControls.H"
#include "CourantNo.H"
//#include "setDeltaT.H"
            runTime.setEndTime(finalTime);
            runTime++;

            Info << "Time = " << runTime.timeName() << nl << endl;

            while (pimple.loop())
            {
                if (pimple.firstIter() || moveMeshOuterCorrectors)
                {
                    // Do any mesh changes
                    mesh.controlledUpdate();
                    std::cerr << "/////////////////////////////////////////////////"<< std::endl;
                 //    sDRBMS.solve();
	                // mesh.movePoints(sDRBMS.curPoints());
	                // std::cerr << "/////////////////////////////////////////////////"<< std::endl;
	                // centerofmassx.append(sDRBMS.motion().centreOfMass().x());
	                // centerofmassy.append(sDRBMS.motion().centreOfMass().y());
	                // centerofmassz.append(sDRBMS.motion().centreOfMass().z());
	                // // To append the linear velocities
	                // velx.append(sDRBMS.motion().v().x());
	                // vely.append(sDRBMS.motion().v().y());
	                // velz.append(sDRBMS.motion().v().z());

                    if (mesh.changing())
                    {
                        MRF.update();

                        if (correctPhi)
                        {
                            // Calculate absolute flux
                            // from the mapped surface velocity
                            phi = mesh.Sf() & Uf();

#include "correctPhi.H"

                            // Make the flux relative to the mesh motion
                            fvc::makeRelative(phi, U);
                        }

                        if (checkMeshCourantNo)
                        {
#include "meshCourantNo.H"
                        }
                    }
                }

#include "UEqn.H"

                // // Solve the Momentum equation
                // MRF.correctBoundaryVelocity(U);

                // tmp<fvVectorMatrix> tUEqn
                // (
                //     fvm::ddt(U) + fvm::div(phi, U)
                //     + MRF.DDt(U)
                //     + turbulence->divDevReff(U)
                //     ==
                //     fvOptions(U)
                // );
                // fvVectorMatrix& UEqn = tUEqn.ref();

                // UEqn.relax();
                // fvOptions.constrain(UEqn);
                // // #############Galerkin projection for the velocity ###########################
                // List<Eigen::MatrixXd> RedLinSysU;
                // if (pimple.momentumPredictor())
                // {
                //     std::cerr << "File: 22fsi.C, Line: 247"<< std::endl;
                //     //solve(UEqn == -fvc::grad(p));
                //     RedLinSysU = Umodes.project(UEqn, NmodesUproj);
                //     volVectorField gradpfull = -fvc::grad(p);
                //     Eigen::MatrixXd projGrad = Umodes.project(gradpfull, NmodesUproj);
                //     RedLinSysU[1] = RedLinSysU[1] + projGrad;
                //     a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual);
                //     Umodes.reconstruct(U, a, "U");
                //     //U.correctBoundaryConditions();
                //     fvOptions.correct(U);
                // }

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    
#include "pEqn.H"

//                     volScalarField rAU(1.0 / UEqn.A());
//                     volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p)); //p
//                     surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
//                     // if (pimple.ddtCorr())
//                     // {
//                     //     phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU) * fvc::ddtCorr(U, phi, Uf));
//                     // }
//                     // else
//                     // {
//                     //     phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU));
//                     // }

//                     //MRF.makeRelative(phiHbyA);

//                     //if (p.needReference())
//                     //{
//                         fvc::makeRelative(phiHbyA, U);
//                         adjustPhi(phiHbyA, U, p);
//                         fvc::makeAbsolute(phiHbyA, U);
//                     //

//                     tmp<volScalarField> rAtU(rAU);

//                     if (pimple.consistent())
//                     {
//                         rAtU = 1.0 / max(1.0 / rAU - UEqn.H1(), 0.1 / rAU);
//                         phiHbyA +=
//                             fvc::interpolate(rAtU() - rAU) * fvc::snGrad(p) * mesh.magSf(); // p
//                         HbyA -= (rAU - rAtU()) * fvc::grad(p); //p
//                     }

//                     if (pimple.nCorrPISO() <= 1)
//                     {
//                         tUEqn.clear();
//                     }

//                     // Update the pressure BCs to ensure flux consistency
//                     // constrainPressure(p, U, phiHbyA, rAtU(), MRF); //p

//                     // ### Reduced linear system for Pressure
//                     List<Eigen::MatrixXd> RedLinSysP;
//                     bold = b;

//                     // Non-orthogonal pressure corrector loop
//                     while (pimple.correctNonOrthogonal())
//                     {
//                         fvScalarMatrix pEqn
//                         (
//                             fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA) //p
//                         );

//                         //pEqn.setReference(pRefCell, pRefValue);

//                         //pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter()))); //p
//                         RedLinSysP = Pmodes.project(pEqn, NmodesPproj);
//                         b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
//                         Pmodes.reconstruct(p, b, "p");

//                         if (pimple.finalNonOrthogonalIter())
//                         {
//                             phi = phiHbyA - pEqn.flux();
//                         }
//                     }

// //#include "continuityErrs.H"

//                     // Explicitly relax pressure for momentum corrector
//                     //p.relax();
//                     //b = bold + mesh.fieldRelaxationFactor("p")*(b-bold);
//                     //Pmodes.reconstruct(p, b, "p");
//                     U = HbyA - rAtU * fvc::grad(p); //p
//                     U.correctBoundaryConditions();
//                     fvOptions.correct(U);
//                     // Correct Uf if the mesh is moving
//                     fvc::correctUf(Uf, U, phi);

//                     // Make the fluxes relative to the mesh motion
//                     fvc::makeRelative(phi, U);

                }// end of the pimple.correct()
                if (pimple.turbCorr())
	            {
	                laminarTransport.correct();
	                turbulence->correct();
	            }

            }// end of the pimple.loop()
                 // the following lines increase the errors
                //Umodes.reconstruct(U, a, "U"); 
                //Pmodes.reconstruct(p, b, "p");
                if(checkWrite(runTime))
                {

                    ITHACAstream::exportSolution(U, name(counter), folder);
                    ITHACAstream::exportSolution(p, name(counter), folder);
                    ITHACAstream::exportSolution(sDRBMS.pointDisplacement(), name(counter), folder);
                    ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
                    std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
                    UredFields.append(U.clone());
			        PredFields.append(p.clone());
			        DFields.append(sDRBMS.pointDisplacement().clone());
			        counter++;
                    nextWrite += writeEvery;
               }
        } // end of the runTime.run() loop

    } // end of the method SolveOnlinePimple

    bool checkWrite(Time& timeObject)
    {
        scalar diffnow = mag(nextWrite - atof(timeObject.timeName().c_str()));
        scalar diffnext = mag(nextWrite - atof(timeObject.timeName().c_str()) -
                              timeObject.deltaTValue());

        if ( diffnow < diffnext)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

};



/*----------------------------------------------------------------------------------------------------------*\
                               Starting the MAIN
\*-----------------------------------------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial22 example(argc, argv);
    tutorial22 online(argc,argv);
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example.meshPtr(),example._runTime());

    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 10);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 10);

    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 5);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 5);

    // Read the par file where the parameters are stored
    word filename("./par");
    example.mu = ITHACAstream::readMatrix(filename);
    // // example.Pnumber = 1; // Number of parameters.
    // // example.Tnumber = 1; //Dimension of the training set (used only when gerating parameters without input)
    // // example.setParameters();
    // // // Set the parameter ranges: Range of the parameter spaces.
    // // example.mu_range(0, 0) = 0.005;
    // // example.mu_range(0, 1) = 0.005;
    // //   // Generate equispaced samples inside the parameter range
    // // example.genEquiPar();

    //Set the inlet boundaries patch 0 directions x and y
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    // Time parameters: We can use Ioodictionnary to access time parameters
    example.startTime = 0;
    example.finalTime = 0.03;
    example.timeStep = 0.01; //0.01;
    example.writeEvery = 0.01;

    // //Perform the offline solve
    example.offlineSolve();
    // //Search the lift function
    // //example.liftSolve3();
    // // Normalize the lifting function
    // //ITHACAutilities::inormalizeFields(example.liftfield);
    // //Create homogeneous basis functions for velocity
    // //example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    // // Perform a POD decomposition for velocity and pressure
    // // ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
    // //                         example.podex, 0, 0, NmodesUout);

    //Perform POD on velocity pressure store the first 20 modes

    // ITHACAPOD::getModes(example.Ufield, online.Umodes, online._U().name(),
    //                     example.podex, 0, 0, NmodesUout);


    // ITHACAPOD::getModes(example.Pfield, online.Pmodes, example._p().name(),
    //                     example.podex, 0, 0,
    //                     NmodesPout);

   //  ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(),
   //                      example.podex, 0, 0, NmodesUout);


   //  ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
   //                      example.podex, 0, 0, NmodesPout);

   //  // Eigen::MatrixXd coeffU = ITHACAutilities::getCoeffs(example.Ufield, online.Umodes, NmodesUproj, true);
   //  // std::cout << "cols of coeffU is :" << coeffU.rows() << std::endl;
    
   //  // Eigen::MatrixXd coeffp= ITHACAutilities::getCoeffs(example.Pfield, online.Pmodes, NmodesPproj, true);
   //  // std::cout << "rows of coeffp is :" << coeffp.rows() << std::endl;

   //  Eigen::MatrixXd coeffU = ITHACAutilities::getCoeffs(example.Ufield, example.Umodes, NmodesUproj, true);
   //  Eigen::MatrixXd coeffp= ITHACAutilities::getCoeffs(example.Pfield, example.Pmodes, NmodesPproj, true);

   // ITHACAstream::exportMatrix(coeffp, "coeffp", "python","./ITHACAoutput/Matrices/");
   // ITHACAstream::exportMatrix(coeffU, "coeffp", "python","./ITHACAoutput/Matrices/");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //example.restart();

    Eigen::VectorXd FoamVely = Foam2Eigen::field2Eigen(example.vely);
    Eigen::VectorXd FoamCentersOfMassy = Foam2Eigen::field2Eigen(example.centerofmassy);

    ITHACAstream::exportMatrix(FoamVely, "FoamVely", "python","./ITHACAoutput/VelOfDispl/");
    ITHACAstream::exportMatrix(FoamCentersOfMassy, "FoamCentersOfMassy", "python","./ITHACAoutput/CenterOfMass/");
    // // ############### contruct the reduced the class object ###################
    //reducedBasicFsi reduced(online);
    reducedBasicFsi reduced(example);
    //example.restart();
    // Reads inlet volocities boundary conditions.
    word vel_file(para->ITHACAdict->lookup("online_velocities"));
    Eigen::MatrixXd vel = ITHACAstream::readMatrix(vel_file);

    reduced.startTime = example.startTime;
    reduced.finalTime = example.finalTime;
    reduced.timeStep = example.timeStep;
    reduced.writeEvery = example.writeEvery;
    //reduced.nextStore = 0.1;
    //reduced.exportEvery = 0.005;

    //Perform the online solutions
    //for (label k = 0; k < (example.mu).size(); k++)
    //{
    //scalar mu_now = example.mu(0, k);
    scalar mu_now = example.mu(0, 0);
    //example.change_viscosity(mu_now);
    //reduced.OnlineVelocity(vel);
    
    //example.restart();
    // ITHACAstream::writePoints(example.meshPtr().points(), "prova", "/polyMesh/");
    reduced.solveOnline_Pimple(mu_now, NmodesUproj, NmodesPproj);
    //}

    // Eigen::MatrixXd errFrobU = ITHACAutilities::errorFrobRel(example.Ufield, reduced.uRecFields);
    // Eigen::MatrixXd errFrobP =  ITHACAutilities::errorFrobRel(example.Pfield,reduced.pRecFields);
    
    // ITHACAstream::exportMatrix(errFrobU, "errFrobU", "matlab","./ITHACAoutput/ErrorsFrob/");
    // ITHACAstream::exportMatrix(errFrobP, "errFrobP", "matlab","./ITHACAoutput/ErrorsFrob/");
   
    // Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Rel(example.Ufield,reduced.uRecFields);
    // Eigen::MatrixXd errL2P =  ITHACAutilities::errorL2Rel(example.Pfield,reduced.pRecFields);
    
    // ITHACAstream::exportMatrix(errL2U, "errL2U", "matlab","./ITHACAoutput/ErrorsL2/");
    // ITHACAstream::exportMatrix(errL2P, "errL2P", "matlab","./ITHACAoutput/ErrorsL2/");
    //scalar mu_now = example.mu(0, 0);
    //example.change_viscosity(mu_now);
    //reduced.OnlineVelocity(vel);
    //reduced.solveOnline_Pimple(mu_now, NmodesUproj, NmodesPproj);
    std::cout << "======================= ONLINE PHASE COMPLETED ================================" << "\n";
    /// Convert Foam List to Eigen Matrices
    Eigen::VectorXd RedVely = Foam2Eigen::field2Eigen(reduced.vely);
    Eigen::VectorXd RedCentersOfMassy = Foam2Eigen::field2Eigen(reduced.centerofmassy);

    //reduced.reconstructLiftAndDrag(CoeffU, Coeffp, "./ITHACAoutput/forces");
    /// exporting matrices
    ITHACAstream::exportMatrix(RedVely, "RedVely", "python","./ITHACAoutput/VelOfDispl/");
    ITHACAstream::exportMatrix(RedCentersOfMassy, "RedCentersOfMassy", "python","./ITHACAoutput/CenterOfMass/");
    
    std::cout << "======================= errorL2Rel ================================" << "\n";
    Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Rel(example.Ufield, reduced.UredFields);
    std::cout << "======================= errL2U completed================================" << "\n";
    Eigen::MatrixXd errL2P = ITHACAutilities::errorL2Rel(example.Pfield, reduced.PredFields);
    std::cout << "======================= errL2P completed================================" << "\n";

    ITHACAstream::exportMatrix(errL2U, "errL2U", "python","./ITHACAoutput/ErrorsL2/");
    ITHACAstream::exportMatrix(errL2P, "errL2P", "python","./ITHACAoutput/ErrorsL2/");

    std::cout << "======================= errorFrobRel ================================" << "\n";
    Eigen::MatrixXd errFrobU = ITHACAutilities::errorFrobRel(example.Ufield, reduced.UredFields);
    std::cout << "======================= errFobU completed================================" << "\n";
    ITHACAstream::exportMatrix(errFrobU, "errFrobU", "python","./ITHACAoutput/ErrorsFrob/");
    Eigen::MatrixXd errFrobP = ITHACAutilities::errorFrobRel(example.Pfield, reduced.PredFields);
    std::cout << "======================= errFobP completed================================" << "\n";

    exit(0);
}


//////////////////////////////////////////////////////////////////////
