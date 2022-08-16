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
    Example of a fluid structure interaction Reduction Problem
SourceFiles
    22fsi.C
\*---------------------------------------------------------------------------*/

#include "fsiBasic.H"
#include "forces.H"
#include "ITHACAPOD.H"
#include "ReducedSimpleSteadyNS.H"
#include "ITHACAstream.H"
#include "dynamicFvMesh.H"
#include "ReducedProblem.H"
#include "Fstream.H"
#include "Foam2Eigen.H"
#include "ITHACAforces.H"
#include "pointVolInterpolation.H"
#include <chrono>
#include<math.h>
#include<iomanip>



class tutorial22: public fsiBasic
{

public:
    explicit tutorial22(int argc, char* argv[])
        : fsiBasic(argc, argv), U(_U()), p(_p()), phi(_phi()), pointDisplacement(_pointDisplacement())
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
    /// PointDisplacement field
    pointVectorField& pointDisplacement;
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
       
         // a = ITHACAutilities::getCoeffs(U, problem->Umodes, 5, true);
         // // std::cout << "cols of a is :" << a.rows() << std::endl;
    
         // b = ITHACAutilities::getCoeffs(p, problem->Pmodes, 5, true);
         // // std::cout << "rows of p is :" << b.rows() << std::endl;
         //uu = problem->Umodes;
         //pp = problem->Pmodes;
        for (int i = 0; i < problem->Umodes.size(); i++)
        {
            uu.append((problem->Umodes.toPtrList()[i]).clone());
        }

        for (int i = 0; i < problem->Pmodes.size(); i++)
        {
            pp.append((problem->Pmodes.toPtrList()[i]).clone());
        }
          //problem->restart();

        std::cout << "################ ctor of reducedBasicFsi ##################" << std::endl;
    }

    tutorial22* problem;

    // Eigen::MatrixXd a;
    // Eigen::MatrixXd b;
    PtrList<volScalarField> PredFields;
    PtrList<volVectorField> UredFields;
    PtrList<pointVectorField> DredFields;
    List<Eigen::MatrixXd> CoeffU;
    List<Eigen::MatrixXd> Coeffp;
    /// Reduced data of centers of mass and velocities displacement
    List<scalar> CenterOfMassx;
    List<scalar> CenterOfMassy;
    List<scalar> CenterOfMassz;
        ///- List of vector: Velocity
    List<scalar> Velx;
    List<scalar> Vely;
    List<scalar> Velz;

    /// Lifted velocity modes.
    volVectorModes uu;
    volScalarModes pp;

    
    ///////// time control variables
    scalar startTime;
    scalar finalTime;
    scalar timeStep;
    scalar writeEvery = timeStep;
    scalar nextWrite;

    label counter = 1;

    /// void solveOnline_Pimple(scalar mu_now, Eigen::MatrixXd& a, Eigen::MatrixXd& b,fileName folder = "./ITHACAoutput/Reconstruct/")
    void solveOnline_Pimple(scalar mu_now, int NmodesUproj, int NmodesPproj,fileName folder = "./ITHACAoutput/Online/")
    {

        // for (int i = 0; i < problem->inletIndex.rows(); i++)
        // {
        //     ULmodes.append((problem->liftfield[i]).clone());
        // }

        // for (int i = 0; i < NmodesUproj; i++)
        // {
        //     ULmodes.append((problem->Umodes.toPtrList()[i]).clone());
        // }

        Eigen::VectorXd uresidual = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::VectorXd presidual = Eigen::VectorXd::Zero(NmodesPproj);

        Eigen::VectorXd a = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::VectorXd b = Eigen::VectorXd::Zero(NmodesPproj);
        problem->restart();
        volScalarField& p = problem->_p();
        volVectorField& U = problem->_U();
        a = ITHACAutilities::getCoeffs(U, uu, NmodesUproj, true);
        b = ITHACAutilities::getCoeffs(p, pp, NmodesPproj, true);
        
        Time& runTime = problem->_runTime();
        surfaceScalarField& phi = problem->_phi();
        dynamicFvMesh& mesh = problem->meshPtr();
        fv::options& fvOptions = problem->_fvOptions();
        pimpleControl& pimple = problem->_pimple();
        pointVectorField& pointDisplacement = problem->_pointDisplacement();
        IOMRFZoneList& MRF = problem->_MRF();
        singlePhaseTransportModel& laminarTransport = problem->_laminarTransport();
        autoPtr<incompressible::turbulenceModel> turbulence =problem->turbulence;
    
        instantList Times = runTime.times();
        runTime.setEndTime(finalTime);
        runTime.setTime(Times[1], 1);
        runTime.setDeltaT(timeStep);
        nextWrite = startTime;
        label pRefCell = 0;
        scalar pRefValue = 0.0;
        scalar  cumulativeContErr =0.0;
        scalar VelDisplX;
        scalar CentreOfMassX;
        IOdictionary dynamicMeshDict
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

        /// Construction of a sixDofRigidBodyMotionSolver Object
        sixDoFRigidBodyMotionSolver sDRBMS(mesh, dynamicMeshDict);
        pointDisplacement = sDRBMS.pointDisplacement();
        const pointMesh& pMesh = pointMesh::New(mesh);
        pointVolInterpolation pvf(pMesh, mesh);
        // Current and old point displacements
        // pointVectorField& displacement(sDRBMS.pointDisplacement());
        // const vectorField displacementOld(mesh().points() - sDRBMS.points0());
        //volVectorField tvf = pvf.interpolate(sDRBMS.pointDisplacement());
   
#include "addCheckCaseOptions.H"
#include "createDyMControls.H"
#include "createUfIfPresent.H"
#include "CourantNo.H"
#include "setInitialDeltaT.H"

        //turbulence->validate();
        //uu.reconstruct(U, a, "U");
        //pp.reconstruct(p, b, "p");
        ITHACAstream::exportSolution(U, name(counter), folder);
        ITHACAstream::exportSolution(p, name(counter), folder);
        ITHACAstream::exportSolution(sDRBMS.pointDisplacement(), name(counter), folder);
        ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
        std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
        UredFields.append(U.clone());
        PredFields.append(p.clone());
        DredFields.append(sDRBMS.pointDisplacement().clone());

        counter++;
        nextWrite += writeEvery;

        // PIMPLE algorithm starts here
        Info<< "\nStarting time loop\n" << endl;
        while (runTime.run())
        {
            std::cerr << "File: 22fsi.C, Line: 230"<< std::endl;
#include "readDyMControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
            runTime.setEndTime(finalTime);
            runTime++;

            Info << "Time = " << runTime.timeName() << nl << endl;
            while (pimple.loop())
            {
                
                if (pimple.firstIter() || moveMeshOuterCorrectors)
                {

                    // std::cerr << "#################### before mesh update ################"<< std::endl;
                    // // // Do any mesh changes
                    //mesh.controlledUpdate();
                    // // //mesh.update();
                    // std::cerr << "################### after mesh update #################"<< std::endl;
                    // std::cerr << "#################### before mesh motionSolver ################"<< std::endl;
                    //mesh.movePoints(sDRBMS.sixDoFRigidBodyMotionSolver::newPoints()); 
                    //volVectorField* Uptr = mesh.getObjectPtr<volVectorField>("U");
                    //U.correctBoundaryConditions();
                    //p.correctBoundaryConditions();
                    // // std::cerr << "################### after mesh motionSolver #################"<< std::endl;
                    /// Solve the Rigid Body motion
                    sDRBMS.solve(); // update the new point displacement
                    // Current point displacements
                    pointVectorField& displacement(sDRBMS.pointDisplacement());
                    // Move the new current points
                    mesh.movePoints(sDRBMS.curPoints());
                    // Old point displacements
                    const vectorField displacementOld(mesh.points() - sDRBMS.points0());
                    vectorField displacementField = displacement.primitiveField() - displacementOld;
                    //Info << displacementField;
                    //ITHACAstream::exportMatrix(displacementField, "displacement"+name(counter), "python", "./ITHACAoutput/displacement/");
                   
                    // if(mesh.getObjectPtr<volVectorField>("U")){
                    //     U.correctBoundaryConditions();
                    // }
                    //U.correctBoundaryConditions();
                    CenterOfMassx.append(sDRBMS.motion().centreOfMass().x());
                    CenterOfMassy.append(sDRBMS.motion().centreOfMass().y());
                    CenterOfMassz.append(sDRBMS.motion().centreOfMass().z()); 
                    ///- List of vector: Velocity
                    //std::cout << "/////////////////"<< sDRBMS.motion().v().y() << std::endl;
                    Velx.append(sDRBMS.motion().v().x());
                    Vely.append(sDRBMS.motion().v().y());
                    Velz.append(sDRBMS.motion().v().z());
                    /// Divide all the velocities by time to have the displacement
                    std::cout << "================= End of the Rigid Body motion Part ============================="<< "\n";
                   
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
                     std::cerr << "File: 22fsi.C, Line: 263"<< std::endl;
                    
                }

                //#include "UEqn.H"

                MRF.correctBoundaryVelocity(U);

                tmp<fvVectorMatrix> tUEqn
                (
                    fvm::ddt(U) + fvm::div(phi, U)
                    + MRF.DDt(U)
                    + turbulence->divDevReff(U)
                    // ==
                    // fvOptions(U)

                );
                fvVectorMatrix& UEqn = tUEqn.ref();

                //UEqn.relax(); // when commenting this line the error drops considerably both U and p

                fvOptions.constrain(UEqn);
                // #############Galerkin projection for the velocity ###########################
                List<Eigen::MatrixXd> RedLinSysU;
                if (pimple.momentumPredictor())
                {
                    std::cerr << "File: 22fsi.C, Line: 247"<< std::endl;
                    //solve(UEqn == -fvc::grad(p));
                    RedLinSysU = uu.project(UEqn, NmodesUproj);
                    volVectorField gradpfull = -fvc::grad(p);
                    Eigen::MatrixXd projGrad = uu.project(gradpfull, NmodesUproj);
                    RedLinSysU[1] = RedLinSysU[1] + projGrad;
                    a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual);
                    uu.reconstruct(U, a, "U");
                    //fvOptions.correct(U);
                    //U.correctBoundaryConditions(); //to check the effect on v()
                }

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    //#include "pEqn.H"

                    volScalarField rAU(1.0 / UEqn.A());
                    volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p)); //p
                    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
                    if (pimple.ddtCorr())
                    {
                        phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU) * fvc::ddtCorr(U, phi, Uf));
                    }
                    else
                    {
                        phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU));
                    }

                    MRF.makeRelative(phiHbyA);

                    if (p.needReference())
                    {
                        fvc::makeRelative(phiHbyA, U);
                        adjustPhi(phiHbyA, U, p);
                        fvc::makeAbsolute(phiHbyA, U);
                    }

                    tmp<volScalarField> rAtU(rAU);

                    if (pimple.consistent())
                    {
                        rAtU = 1.0 / max(1.0 / rAU - UEqn.H1(), 0.1 / rAU);
                        phiHbyA +=
                            fvc::interpolate(rAtU() - rAU) * fvc::snGrad(p) * mesh.magSf(); // p
                        HbyA -= (rAU - rAtU()) * fvc::grad(p); //p
                    }

                    if (pimple.nCorrPISO() <= 1)
                    {
                        tUEqn.clear();
                    }

                    // Update the pressure BCs to ensure flux consistency
                    constrainPressure(p, U, phiHbyA, rAtU(), MRF); //p

                    // ### Reduced linear system for Pressure
                    List<Eigen::MatrixXd> RedLinSysP;
                    //bOld = b;
                    // Non-orthogonal pressure corrector loop
                    while (pimple.correctNonOrthogonal())
                    {
                        fvScalarMatrix pEqn
                        (
                            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA) //p
                        );

                        pEqn.setReference(pRefCell, pRefValue);

                        //pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter()))); //p
                        RedLinSysP = pp.project(pEqn, NmodesPproj);
                        b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
                        pp.reconstruct(p, b, "p");
                        if (pimple.finalNonOrthogonalIter())
                        {
                            phi = phiHbyA - pEqn.flux();
                        }
                    }
                    
                    // Explicitly relax pressure for momentum corrector
                    //p.relax();  // no effect on the reduced field
                    U = HbyA - rAtU * fvc::grad(p); //p
                    U.correctBoundaryConditions();
                    fvOptions.correct(U);
                    // Correct Uf if the mesh is moving
                    fvc::correctUf(Uf, U, phi);

                    // Make the fluxes relative to the mesh motion
                    fvc::makeRelative(phi, U);

                }// out of the pimple.correct() loop

            }// out of the pimple.loop() loop
            uu.reconstruct(U, a, "U");
            pp.reconstruct(p, b, "p");
            
            if(checkWrite(runTime))
            {
                
                ITHACAstream::exportSolution(U, name(counter), folder);
                ITHACAstream::exportSolution(p, name(counter), folder);
                ITHACAstream::exportSolution(sDRBMS.pointDisplacement(), name(counter), folder);
                ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
              
                std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
                UredFields.append(U.clone());
                PredFields.append(p.clone());
                DredFields.append(sDRBMS.pointDisplacement().clone());

                //ITHACAstream::exportMatrix(Foam2Eigen::field2Eigen(sDRBMS.pointDisplacement()), "pd"+counter, "python","./ITHACAoutput/pythonDispl/");
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

/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/

int main(int argc, char* argv[])
{
    // Construct the tutorial object
    tutorial22 example(argc, argv);

    tutorial22 online(argc,argv);
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example.meshPtr(),example._runTime());

    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 50);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 50);

    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);

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
    example.finalTime = 0.5;
    example.timeStep = 0.005; //0.01;
    example.writeEvery = 0.005;

    // //Perform the offline solve
    example.offlineSolve();

    // //Search the lift function
    // //example.liftSolve3();
    // // Normalize the lifting function
    // //ITHACAutilities::inormalizeFields(example.liftfield);
    // //Create homogeneous basis functions for velocity
    // //example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    // // Perform a POD decomposition for velocity and pressure including lifting function
    // // ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
    // //                         example.podex, 0, 0, NmodesUout);

    //Perform POD on velocity pressure store the first 20 modes

    //ITHACAPOD::getModes(example.Ufield, online.Umodes, online._U().name(),example.podex, 0, 0, NmodesUout);

    //ITHACAPOD::getModes(example.Pfield, online.Pmodes, example._p().name(),example.podex, 0, 0, NmodesPout);
    

    ITHACAPOD::getModes(example.Ufield, example.Umodes, example._U().name(), example.podex, 0, 0, NmodesUout);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(), example.podex, 0, 0, NmodesPout);

    //example.pmodes = example.Pmodes;
    //example.umodes = example.Umodes;

    // ITHACAPOD::getModes(example.Pfield, online.Pmodes, example._p().name(),
    //                     example.podex, 0, 0, NmodesPout);
     // ITHACAPOD::getModes(example.PointDispfield, online.PDmodes, example._pointDisplacement().name(),
     //                    example.podex, 0, 0,
     //                    NmodesPout);


    // Eigen::MatrixXd CoeffU = ITHACAutilities::getCoeffs(example.Ufield, online.Umodes, NmodesUproj, true);
    // //std::cout << "the number rows of U is :" << CoeffU.rows() << "========"<< "the number cols of U is :"<< CoeffU.cols() << "\n";
    // ITHACAstream::exportMatrix  (CoeffU,"coeffUMatrix", "python","./Matrices");

    
    // Eigen::MatrixXd Coeffp= ITHACAutilities::getCoeffs(example.Pfield, online.Pmodes, NmodesPproj, true);
   
    //ITHACAstream::exportMatrix  (Coeffp, "coeffpMatrix", "python","./Matrices");
    //std::cout << "the number rows of p is :" << Coeffp.rows() << "========"<< "the number cols of p is :"<< Coeffp.cols() << "\n";
    Eigen::VectorXd FoamVely = Foam2Eigen::field2Eigen(example.Vely);
    ITHACAstream::exportMatrix(FoamVely, "FoamVely", "python","./ITHACAoutput/VelOfDispl/");
    //Eigen::VectorXd FoamVely = Foam2Eigen::field2Eigen(online.Vely);
    Eigen::VectorXd FoamCentersOfMassy = Foam2Eigen::field2Eigen(example.CenterOfMassy);
    ITHACAstream::exportMatrix(FoamCentersOfMassy, "FoamCentersOfMassy", "python","./ITHACAoutput/CenterOfMass/");
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //example.restart();
    // ############### contruct the reduced the class object ###################
    reducedBasicFsi reduced(example);
 
    // Reads inlet volocities boundary conditions.
    word vel_file(para->ITHACAdict->lookup("online_velocities"));
    Eigen::MatrixXd vel = ITHACAstream::readMatrix(vel_file);

    reduced.startTime =  0;
    reduced.finalTime =  0.5;
    reduced.timeStep =   0.005;
    reduced.writeEvery = 0.005;

    //Perform the online solutions
    //for (label k = 0; k < (example.mu).size(); k++)
    //{
    //scalar mu_now = example.mu(0, k);
    scalar mu_now = example.mu(0, 0);
    //example.change_viscosity(mu_now);
    //reduced.OnlineVelocity(vel);
    reduced.solveOnline_Pimple(mu_now, NmodesUproj, NmodesPproj);
    std::cout << "======================= ONLINE PHASE COMPLETED ================================" << "\n";
    /// Convert Foam List to Eigen Matrices
    Eigen::VectorXd RedVely = Foam2Eigen::field2Eigen(reduced.Vely);
    Eigen::VectorXd RedCentersOfMassy = Foam2Eigen::field2Eigen(reduced.CenterOfMassy);

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
    ITHACAstream::exportMatrix(errFrobP, "errFrobP", "python","./ITHACAoutput/ErrorsFrob/");
    exit(0);
}






// void forcesMatrices(label nModes)
// {
//     tauMatrix.resize(nModes, 3);
//     nMatrix.resize(nModes, 3);
//     tauMatrix = tauMatrix * 0;
//     nMatrix = nMatrix * 0;
//     //Time& runTime = _runTime();
//     Time& runTime = problem->_runTime();
//     instantList Times = runTime.times();
//     //fvMesh& mesh = _mesh();
//     dynamicFvMesh& mesh = problem->meshPtr();
//     surfaceScalarField& phi = problem->_phi();
//     volScalarField& p = problem->_p();
//     volVectorField& U = problem->_U();
//     //Read FORCESdict
//     IOdictionary FORCESdict
//     (
//         IOobject
//         (
//             "FORCESdict",
//             runTime.system(),
//             mesh,
//             IOobject::MUST_READ,
//             IOobject::NO_WRITE
//         )
//     );
//     IOdictionary transportProperties
//     (
//         IOobject
//         (
//             "transportProperties",
//             runTime.constant(),
//             mesh,
//             IOobject::MUST_READ,
//             IOobject::NO_WRITE
//         )
//     );
//     word pName(FORCESdict.lookup("pName"));
//     word UName(FORCESdict.lookup("UName"));
//     functionObjects::ITHACAforces f("Forces", mesh, FORCESdict);
 
//     for (label i = 0; i < nModes; i++)
//     {
//         U = problem->Umodes[i];
//         p = problem->Pmodes[0];
//         mesh.readUpdate();
//         f.write();
//         f.calcForcesMoment();
 
//         for (label j = 0; j < 3; j++)
//         {
//             tauMatrix(i, j) = f.forceTau()[j];
//         }
//     }
 
//     for (label i = 0; i < nModes; i++)
//     {
//         U = problem->Umodes[0];
//         p = problem->Pmodes[i];
//         mesh.readUpdate();
//         f.write();
//         f.calcForcesMoment();
 
//         for (label j = 0; j < 3; j++)
//         {
//             nMatrix(i, j) = f.forcePressure()[j];
//         }
//     }
 
//     if (problem->para->exportPython)
//     {
//         ITHACAstream::exportMatrix(tauMatrix, "tau", "python",
//                                    "./ITHACAoutput/Matrices/");
//         ITHACAstream::exportMatrix(nMatrix, "n", "python", "./ITHACAoutput/Matrices/");
//     }
 
//     if (problem->para->exportMatlab)
//     {
//         ITHACAstream::exportMatrix(tauMatrix, "tau", "matlab",
//                                    "./ITHACAoutput/Matrices/");
//         ITHACAstream::exportMatrix(nMatrix, "n", "matlab", "./ITHACAoutput/Matrices/");
//     }
 
//     if (problem->para->exportTxt)
//     {
//         ITHACAstream::exportMatrix(tauMatrix, "tau", "eigen",
//                                    "./ITHACAoutput/Matrices/");
//         ITHACAstream::exportMatrix(nMatrix, "n", "eigen", "./ITHACAoutput/Matrices/");
//     }

// }


// void reconstructLiftAndDrag(const Eigen::MatrixXd& velCoeffs,
//                                       const Eigen::MatrixXd& pressureCoeffs, fileName folder)
// {
//     M_Assert(velCoeffs.cols() == tauMatrix.rows(),
//              "The number of velocity modes in the coefficients matrix is not equal to the number of modes in the viscous forces matrix.");
//     M_Assert(pressureCoeffs.cols() == nMatrix.rows(),
//              "The number of pressure modes in the coefficients matrix is not equal to the number of modes in the pressure forces matrix.");
//     mkDir(folder);
//     system("ln -s ../../constant " + folder + "/constant");
//     system("ln -s ../../0 " + folder + "/0");
//     system("ln -s ../../system " + folder + "/system");
//     //Read FORCESdict
//     IOdictionary FORCESdict
//     (
//         IOobject
//         (
//             "FORCESdict",
//             "./system",
//             problem->Umodes[0].mesh(),
//             IOobject::MUST_READ,
//             IOobject::NO_WRITE
//         )
//     );
//     Eigen::MatrixXd fTau;
//     Eigen::MatrixXd fN;
//     fTau.setZero(velCoeffs.rows(), 3);
//     fN.setZero(pressureCoeffs.rows(), 3);
//     fTau = velCoeffs * tauMatrix;
//     fN = pressureCoeffs * nMatrix;
 
//     // Export the matrices
//     if (problem->para->exportPython)
//     {
//         ITHACAstream::exportMatrix(fTau, "fTau", "python", folder);
//         ITHACAstream::exportMatrix(fN, "fN", "python", folder);
//     }
 
//     if (problem->para->exportMatlab)
//     {
//         ITHACAstream::exportMatrix(fTau, "fTau", "matlab", folder);
//         ITHACAstream::exportMatrix(fN, "fN", "matlab", folder);
//     }
 
//     if (problem->para->exportTxt)
//     {
//         ITHACAstream::exportMatrix(fTau, "fTau", "eigen", folder);
//         ITHACAstream::exportMatrix(fN, "fN", "eigen", folder);
//     }
// }
