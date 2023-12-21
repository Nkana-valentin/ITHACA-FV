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
#include "dynamicMotionSolverFvMesh.H"
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
        point0 = meshPtr().points();
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

        truthSolve(mu_now);
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
       
        for (int i = 0; i < problem->inletIndex.rows(); i++)
        {
             uu.append((problem->liftfield[i]).clone());
        }
        for (int i = 0; i < problem->Umodes.size(); i++)
        {
            uu.append((problem->Umodes.toPtrList()[i]).clone());
        }

        for (int i = 0; i < problem->Pmodes.size(); i++)
        {
            pp.append((problem->Pmodes.toPtrList()[i]).clone());
        }
        //problem->restart();
    }

    tutorial22* problem;
    PtrList<volScalarField> PredFields;
    PtrList<volVectorField> UredFields;
    PtrList<surfaceScalarField> Phiredfield;
    PtrList<pointVectorField> DredFields;
    PtrList< volVectorField > Ufield;
    PtrList<volScalarField>   Pfield;
    PtrList<pointVectorField> Dfields;
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
    List<scalar> romforcey;
    List<scalar> romforcex;
    ///////// time control variables
    scalar startTime;
    scalar finalTime;
    scalar timeStep;
    scalar writeEvery = timeStep;
    scalar nextWrite;
    label counter = 1;

    void solveOnline_Pimple(scalar mu_now, int NmodesUproj, int NmodesPproj,fileName folder = "./ITHACAoutput/Online/")
    {

        Eigen::VectorXd uresidual = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::VectorXd presidual = Eigen::VectorXd::Zero(NmodesPproj);
        Eigen::VectorXd a = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::MatrixXd a0 = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::VectorXd b = Eigen::VectorXd::Zero(NmodesPproj);
        Eigen::VectorXd bOld = Eigen::VectorXd::Zero(NmodesPproj);
  
        Time& runTime = problem->_runTime();
        dynamicFvMesh& mesh = problem->meshPtr();
        fv::options& fvOptions = problem->_fvOptions();
        pimpleControl& pimple = problem->_pimple();
        volScalarField& p = problem->_p();
        volScalarField pOld = p;
        volVectorField& U = problem->_U();
        pointVectorField& pointDisplacement = problem->_pointDisplacement(); //const_cast<pointVectorField&>(mesh.lookupObject<pointVectorField>("pointDisplacement"));
//ITHACAstream::writePoints(problem->GridPoints[0], folder, name(counter) + "/polyMesh/"); //???
        //ITHACAstream::read_fields(Dfields, pointDisplacement, "./ITHACAoutput/Offline/");
        //std::cout << "==============" << Dfields.size() << "===========" << std::endl;
        ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
        std::cout << "==============" << Ufield.size() << "===========" << std::endl;
        ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
        std::cout << "==============" << Pfield.size() << "===========" << std::endl;
        U = Ufield[0];
        p = Pfield[0];

        sixDoFRigidBodyMotionSolver& sDRBMS = problem->sDRBMS();
        // pointDisplacement.primitiveFieldRef() = problem->Dfield[0].primitiveField();
        // sDRBMS.pointDisplacement().primitiveFieldRef() = pointDisplacement.primitiveFieldRef(); //problem->Dfield[0];

        surfaceScalarField& phi = problem->_phi();
        IOMRFZoneList& MRF = problem->_MRF();
        singlePhaseTransportModel& laminarTransport = problem->_laminarTransport();
        autoPtr<incompressible::turbulenceModel> turbulence = problem->turbulence;
        dictionary dictCoeffs(problem->dyndict->findDict("sixDoFRigidBodyMotionCoeffs"));
        Foam::functionObjects::forces romforces("romforces", mesh,  dictCoeffs);
        turbulence->validate();
#include "createUfIfPresent.H"
        bool  correctPhi = problem->correctPhi;
        bool  checkMeshCourantNo = problem->checkMeshCourantNo;
        bool  moveMeshOuterCorrectors = problem->moveMeshOuterCorrectors;
    
        instantList Times = runTime.times();
        runTime.setEndTime(finalTime);
        runTime.setTime(Times[1], 1);
        runTime.setDeltaT(timeStep);
        nextWrite = startTime;
        label pRefCell = 0;
        scalar pRefValue = 0.0;
        scalar  cumulativeContErr = problem->cumulativeContErr;
        a = ITHACAutilities::getCoeffs(U, uu, NmodesUproj, true);
        b = ITHACAutilities::getCoeffs(p, pp, NmodesPproj, true);
        // Current and old point displacements
        // pointVectorField& displacement(sDRBMS.pointDisplacement());
        // const vectorField displacementOld(mesh().points() - sDRBMS.points0());
        //volVectorField tvf = pvf.interpolate(sDRBMS.pointDisplacement());
        //a(0) = a0(0);
        // Current and old point displacements
        // pointVectorField& displacement(sDRBMS.pointDisplacement());
        // const vectorField displacementOld(mesh().points() - sDRBMS.points0());
        
//#include "createDyMControls.H"
//#include "createUfIfPresent.H"
//#include "CourantNo.H"
//#include "setInitialDeltaT.H"
        uu.reconstruct(U, a, "U");
        pp.reconstruct(p, b, "p");

        // PIMPLE algorithm starts here
        Info<< "\nStarting time loop\n" << endl;
        while (runTime.run())
        {
        
//#include "readDyMControls.H"
//#include "CourantNo.H"
//#include "setDeltaT.H"
            runTime.setEndTime(finalTime);
            runTime++;

            Info << "Time = " << runTime.timeName() << nl << endl;
            while (pimple.loop())
            {
                
                if (pimple.firstIter() || moveMeshOuterCorrectors)
                {
                    //mesh.controlledUpdate();
                    romforces.calcForcesMoment();
                    // // Current point displacements
                    // pointVectorField& displacement(sDRBMS.pointDisplacement());
                    // /// Solve the Rigid Body motion
                    sDRBMS.solve(); // update the new point displacement
                    // // // Move the new current points
                     mesh.movePoints(sDRBMS.curPoints());
                
                    if (mesh.changing())
                    {
                        MRF.update();

                        if (correctPhi)
                        {
                            // Calculate absolute flux
                            // from the mapped surface velocity
                            phi = mesh.Sf() & Uf();

//#include "correctPhi.H"
                            // Make the flux relative to the mesh motion
                            fvc::makeRelative(phi, U);
                        }

                        if (checkMeshCourantNo)
                        {
#include "meshCourantNo.H"
                        }
                    }
                    
                }        

                MRF.correctBoundaryVelocity(U);

                tmp<fvVectorMatrix> tUEqn
                (
                    fvm::ddt(U) 
                    + MRF.DDt(U)
                    + turbulence->divDevReff(U)
                    == -fvc::div(phi, U) + fvOptions(U)

                );
                fvVectorMatrix& UEqn = tUEqn.ref();

                UEqn.relax(); 
                fvOptions.constrain(UEqn);
                List<Eigen::MatrixXd> RedLinSysU;
                if (pimple.momentumPredictor())
                {
                    // solve(UEqn == -fvc::grad(p));
                    RedLinSysU = uu.project(UEqn, NmodesUproj);
                    volVectorField gradpfull = -fvc::grad(p);
                    Eigen::MatrixXd projGrad = uu.project(gradpfull, NmodesUproj);
                    RedLinSysU[1] = RedLinSysU[1] + projGrad;
                    a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual);
                    uu.reconstruct(U, a, "U");
                    fvOptions.correct(U); //?
                }
                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    volScalarField rAU(1.0 / UEqn.A());
                    volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));
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
                    constrainPressure(p, U, phiHbyA, rAtU(), MRF);
                    List<Eigen::MatrixXd> RedLinSysP;
                    bOld = b;
                    // Non-orthogonal pressure corrector loop
                    while (pimple.correctNonOrthogonal())
                    {
                        fvScalarMatrix pEqn
                        (
                            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
                        );

                        pEqn.setReference(pRefCell, pRefValue);
                        RedLinSysP = pp.project(pEqn, NmodesPproj);
                        b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
                        //b = bOld + mesh.fieldRelaxationFactor("p") * (b - bOld); // relax with reduced coeffs
                        pp.reconstruct(p, b, "p");
                        if (pimple.finalNonOrthogonalIter())
                        {
                            phi = phiHbyA - pEqn.flux();
                        }
                    }
                    // relax pressure for momentum corrector
                    //b = bOld + mesh.fieldRelaxationFactor("p") * (b - bOld);
                    //pp.reconstruct(p, b, "p");
                    //#include "continuityErrs.H"

                    // Explicitly relax pressure for momentum corrector
                    p.relax();
                    U = HbyA - rAtU * fvc::grad(p);
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
                romforcey.append(romforces.forceEff().y());
                romforcex.append(romforces.forceEff().x());

                CenterOfMassx.append(sDRBMS.motion().centreOfMass().x());
                CenterOfMassy.append(sDRBMS.motion().centreOfMass().y());
                CenterOfMassz.append(sDRBMS.motion().centreOfMass().z()); 
                ///- List of vector: Velocity
                Velx.append(sDRBMS.motion().v().x());
                Vely.append(sDRBMS.motion().v().y());
                Velz.append(sDRBMS.motion().v().z());  
                
                ITHACAstream::exportSolution(U, name(counter), folder);
                ITHACAstream::exportSolution(p, name(counter), folder);
                ITHACAstream::exportSolution(phi, name(counter), folder);
                ITHACAstream::exportSolution(sDRBMS.pointDisplacement(), name(counter), folder);
                ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
              
                std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
                UredFields.append(U.clone());
                PredFields.append(p.clone());
                Phiredfield.append(phi.clone());
                DredFields.append(sDRBMS.pointDisplacement().clone());
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

    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj",25);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj",25);

    // Read the par file where the parameters are stored
    word filename("./par");
    example.mu = ITHACAstream::readMatrix(filename);

    //Set the inlet boundaries patch 0 directions x and y
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;

    // Time parameters: We can use Ioodictionnary to access time parameters
    example.startTime  = 50;
    example.finalTime  = 100;
    example.timeStep   = 5e-03;
    example.writeEvery = 5e-01;

    // //Perform the offline solve
    example.offlineSolve();

    std::cout << "==============" << example.Dfield.size() << "==============" << std::endl;
    online.Dfield = example.Dfield;
    //Search the lift function
    //example.liftSolve();
    online.liftSolve();
    online.GridPoints = example.GridPoints;

    ITHACAPOD::getModes(example.Ufield, online.Umodes, online._U().name(),example.podex, 0, 0, NmodesUout);

    ITHACAPOD::getModes(example.Pfield, online.Pmodes, example._p().name(),example.podex, 0, 0, NmodesPout);

    Eigen::MatrixXd CoeffU = ITHACAutilities::getCoeffs(example.Ufield, online.Umodes, NmodesUproj, true);
    Eigen::MatrixXd Coeffp = ITHACAutilities::getCoeffs(example.Pfield, online.Pmodes, NmodesPproj, true);
   
    ITHACAstream::exportMatrix  (CoeffU, "coeffUMatrix", "python","./Matrices");
    ITHACAstream::exportMatrix  (Coeffp, "coeffpMatrix", "python","./Matrices");

    Eigen::VectorXd FoamVely = Foam2Eigen::field2Eigen(example.vely);
    ITHACAstream::exportMatrix(FoamVely, "FoamVely", "python","./ITHACAoutput/VelOfDispl/");
    Eigen::VectorXd FoamCentersOfMassy = Foam2Eigen::field2Eigen(example.centerofmassy);
    ITHACAstream::exportMatrix(FoamCentersOfMassy, "FoamCentersOfMassy", "python","./ITHACAoutput/CenterOfMass/");

    Eigen::VectorXd fomforcex = Foam2Eigen::field2Eigen(example.fomforcex);
    ITHACAstream::exportMatrix(fomforcex, "fomforcex", "python","./ITHACAoutput/forces/");

    Eigen::VectorXd fomforcey = Foam2Eigen::field2Eigen(example.fomforcey);
    ITHACAstream::exportMatrix(fomforcey, "fomforcey", "python","./ITHACAoutput/forces/");
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ############### contruct the reduced the class object ###################
    // reducedBasicFsi reduced(example);
     reducedBasicFsi reduced(online);

    reduced.startTime =  example.startTime;
    reduced.finalTime =  example.finalTime;
    reduced.timeStep =   example.timeStep;
    reduced.writeEvery = example.writeEvery;
    scalar mu_now = example.mu(0, 0);
    //example.change_viscosity(mu_now);
    //reduced.OnlineVelocity(vel);
    reduced.solveOnline_Pimple(mu_now, NmodesUproj, NmodesPproj);
    std::cout << "======================= ONLINE PHASE COMPLETED ================================" << "\n";
    /// Convert Foam List to Eigen Matrices
    Eigen::VectorXd RedVely = Foam2Eigen::field2Eigen(reduced.Vely);
    Eigen::VectorXd RedCentersOfMassy = Foam2Eigen::field2Eigen(reduced.CenterOfMassy);

    /// exporting matrices
    ITHACAstream::exportMatrix(RedVely, "RedVely", "python","./ITHACAoutput/VelOfDispl/");
    ITHACAstream::exportMatrix(RedCentersOfMassy, "RedCentersOfMassy", "python","./ITHACAoutput/CenterOfMass/");
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Eigen::VectorXd romforcex = Foam2Eigen::field2Eigen(reduced.romforcex);
    ITHACAstream::exportMatrix(romforcex, "romforcex", "python","./ITHACAoutput/forces/");

    Eigen::VectorXd romforcey = Foam2Eigen::field2Eigen(reduced.romforcey);
    ITHACAstream::exportMatrix(romforcey, "romforcey", "python","./ITHACAoutput/forces/");
    
    std::cout << "======================= errorL2Rel ================================" << "\n";
    Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Rel(example.Ufield, reduced.UredFields);
    std::cout << "======================= errL2U completed================================" << "\n";
    Eigen::MatrixXd errL2P = ITHACAutilities::errorL2Rel(example.Pfield, reduced.PredFields);
    std::cout << "======================= errL2P completed================================" << "\n";
    //Eigen::MatrixXd errL2Phi = ITHACAutilities::errorL2Rel(example.Phifield, reduced.Phiredfield);
    //std::cout << "======================= errL2Phi completed================================" << "\n";

    ITHACAstream::exportMatrix(errL2U, "errL2U", "python","./ITHACAoutput/ErrorsL2/");
    ITHACAstream::exportMatrix(errL2P, "errL2P", "python","./ITHACAoutput/ErrorsL2/");

    std::cout << "======================= errorFrobRel ================================" << "\n";
    Eigen::MatrixXd errFrobU = ITHACAutilities::errorFrobRel(example.Ufield, reduced.UredFields);
    std::cout << "======================= errFobU completed================================" << "\n";
    ITHACAstream::exportMatrix(errFrobU, "errFrobU", "python","./ITHACAoutput/ErrorsFrob/");
    Eigen::MatrixXd errFrobP = ITHACAutilities::errorFrobRel(example.Pfield, reduced.PredFields);
    std::cout << "======================= errFobP completed================================" << "\n";
    //Eigen::MatrixXd errFrobPhi = ITHACAutilities::errorFrobRel(example.Phifield, reduced.Phiredfield);
    //std::cout << "======================= errFobP completed================================" << "\n";
    //ITHACAstream::exportMatrix(errFrobP, "errFrobP", "python","./ITHACAoutput/ErrorsFrob/");
    exit(0);
}
