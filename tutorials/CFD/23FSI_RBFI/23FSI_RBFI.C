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
#include "pointConstraints.H"
#include "zoneMotion.H"

class tutorial22: public fsiBasic
{
public:
    explicit tutorial22(int argc, char* argv[])
        : fsiBasic(argc, argv), U(_U()), p(_p()), phi(_phi()), pd(_pd())
    {
        //point0 = meshPtr().points();
        //restart();
    }

    // Fields To Perform
    /// Velocity field
    volVectorField& U;
    /// Pressure field
    volScalarField& p;
    /// flux field
    surfaceScalarField& phi;
    pointVectorField& pd;
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
            // for (int i = 0; i < problem->inletIndex.rows(); i++)
            // {
            //     Umodes.append((problem->liftfield[i]).clone());
            // }
            for (int i = 0; i < problem->Umodes.size(); i++)
		    {
		        Umodes.append((problem->Umodes.toPtrList()[i]).clone());
		    }

		    for (int i = 0; i < problem->Pmodes.size(); i++)
		    {
		        Pmodes.append((problem->Pmodes.toPtrList()[i]).clone());
		    }

            for (int i = 0; i < problem->Dmodes.size(); i++)
            {
                Dmodes.append((problem->Dmodes.toPtrList()[i]).clone());
            }
		    //problem->restart();

        std::cout << "################ ctor of POD-I Fsi ##################" << std::endl;
    }

    tutorial22* problem;
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

   PtrList<volScalarField> PredFields;
   PtrList<volVectorField> UredFields;
   PtrList<pointVectorField> Dfield;

   label counter = problem->counter;
    // /// Create a samples for interpolation
    // std::vector<SPLINTER::DataTable*> samples;

    // /// Create a RBF splines for interpolation
    // std::vector<SPLINTER::RBFSpline*> rbfSplines;

    void solveOnline_Pimple(scalar mu_now, int NmodesUproj, int NmodesPproj, int NmodesDproj, fileName folder = "./ITHACAoutput/Reconstruct/")
    {
        
    
        for (int i = 0; i < NmodesUproj; i++)
        {
            Umodes.append((problem->Umodes.toPtrList()[i]).clone());
        }

        for (int i = 0; i < NmodesPproj; i++)
        {
            Pmodes.append((problem->Pmodes.toPtrList()[i]).clone());
        }
        Eigen::VectorXd uresidual = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::VectorXd presidual = Eigen::VectorXd::Zero(NmodesPproj);
        //std::cout <<  presidual.cwiseAbs() << std::endl;
        scalar U_norm_res(1);
        scalar P_norm_res(1);

        Eigen::MatrixXd a = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::MatrixXd a0 = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::MatrixXd b = Eigen::VectorXd::Zero(NmodesPproj);
        Eigen::MatrixXd bold = Eigen::VectorXd::Zero(NmodesPproj);
        Eigen::MatrixXd c = Eigen::VectorXd::Zero(NmodesDproj);

        problem->restart();
        Time& runTime = problem->_runTime();
        surfaceScalarField& phi = problem->_phi();
        dynamicFvMesh& mesh = problem->meshPtr();
        fv::options& fvOptions = problem->_fvOptions();
        pimpleControl& pimple = problem->_pimple();
        volScalarField& p = problem->_p();
        volVectorField& U = problem->_U();
        //pointVectorField& pointDisplacement = problem->_pd();//????
        pointVectorField& pointDisplacement = const_cast<pointVectorField&>
                                          (mesh.lookupObject<pointVectorField>("pointDisplacement"));
        //surfaceVectorField& Uf = problem->_Uf();
        IOMRFZoneList& MRF = problem->_MRF();
        singlePhaseTransportModel& laminarTransport = problem->_laminarTransport();
        autoPtr<incompressible::turbulenceModel> turbulence = problem->turbulence;
        instantList Times = runTime.times();
        runTime.setEndTime(finalTime);
        // Declare modal coefficients for velocity and pressure
        a = ITHACAutilities::getCoeffs(U, Umodes, NmodesUproj, true);
        b = ITHACAutilities::getCoeffs(p, Pmodes, NmodesPproj, true);
        c = ITHACAutilities::getCoeffs(pointDisplacement, Dmodes, NmodesDproj, false);
        // To solve the rbfi system to obtain the weights
        //PodIpointDispl(problem->coeffL2, problem->CylDispl, NmodesDproj);
        problem->PodIpointDispl(problem->CylDispl, NmodesDproj);
    
        runTime.setTime(Times[1], 1);
        runTime.setDeltaT(timeStep);
        nextWrite = startTime;
        label pRefCell = 0;
        scalar pRefValue = 0.0;

        bool  correctPhi = problem->correctPhi;
        bool  checkMeshCourantNo = problem->checkMeshCourantNo;
        bool  moveMeshOuterCorrectors = problem->moveMeshOuterCorrectors;
        scalar  cumulativeContErr = problem->cumulativeContErr;
        const scalar ramp = 1.0;
    
    //std::cerr << "/////////////////////"<< problem->sDRBMS().motion().mass() << "///////////////////" << std::endl;
#include "createUfIfPresent.H"
        turbulence->validate();
        //Umodes.reconstruct(U, a, "U");
        //Pmodes.reconstruct(p, b, "p");
        problem->Dmodes.reconstruct(pointDisplacement, c, "pointDisplacement");
        pointDisplacement.primitiveFieldRef() = pointDisplacement;
        ITHACAstream::exportSolution(pointDisplacement, name(counter), folder);
        ITHACAstream::exportSolution(U, name(counter), folder);
        ITHACAstream::exportSolution(p, name(counter), folder);
        ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
        std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
        UredFields.append(U.clone());
        PredFields.append(p.clone());
        counter++;
        nextWrite += writeEvery;
        dictionary dictCoeffs(problem->dyndict->findDict("sixDoFRigidBodyMotionCoeffs"));
        Foam::functionObjects::forces romforces("romforces", mesh, dictCoeffs);
        sixDoFRigidBodyMotion sDRBM(dictCoeffs, dictCoeffs, runTime );
        Foam::dimensionedVector g("g", dimAcceleration, Zero);
        dictCoeffs.readIfPresent("g", g);
        Eigen::MatrixXd pdCoeff;
        pdCoeff.resize(NmodesDproj, 1);
        bool firstIter = false;
        pointField points0 = mesh.points();
        // PIMPLE algorithm starts here
        Info<< "\nStarting time loop\n" << endl;
        while (runTime.run())
        {

            runTime.setEndTime(finalTime);
            runTime++;

            Info << "Time = " << runTime.timeName() << nl << endl;

            while (pimple.loop())
            {
                if (pimple.firstIter() || moveMeshOuterCorrectors)
                {


                    sDRBM.newTime();
                    romforces.calcForcesMoment();
                    romforcex.append(romforces.forceEff().x());
                    romforcey.append(romforces.forceEff().y());
                    // Solving the sixRigidMotion problem
                    sDRBM.update
                    (
                         firstIter,
                         ramp*(romforces.forceEff() + sDRBM.mass()*g.value()),
                         ramp
                        *(
                            romforces.momentEff()
                          + sDRBM.mass()*(sDRBM.momentArm() ^ g.value())
                         ),
                         runTime.deltaTValue(),
                         runTime.deltaT0Value()
                    );

                    //std::cout << "/////////////" <<  sDRBM.centreOfMass().y() << "////////////" << std::endl;
                    for (int i = 0; i < NmodesDproj; i++)
                    {
                        Eigen::MatrixXd muEval;
                        muEval.resize(1, 1);
                        // mu_now = sDRBM.centreOfMass().y();
                        muEval(0, 0)  = sDRBM.centreOfMass().y();
                        pdCoeff(i, 0) = problem->rbfSplines[i]->eval(muEval);
                    }
                    pointConstraints::New(pointDisplacement.mesh()).constrainDisplacement(pointDisplacement);
                    // Reconstruction of the pointdisplacement
                    problem->Dmodes.reconstruct(pointDisplacement, pdCoeff, "pointDisplacement");
                    // Update the displacements
                    //pointDisplacement.primitiveFieldRef() = points0 + pointDisplacement.primitiveField();
                    //pointDisplacement.correctBoundaryConditions();
                    // Displacement has changed. Update boundary conditions
                    //pointConstraints::New(pointDisplacement.mesh()).constrainDisplacement(pointDisplacement);
                    // tmp<pointField> newPoints(mesh.points() + pointDisplacement.primitiveField());
                    tmp<pointField> newPoints(points0 + pointDisplacement.primitiveField());
                    //pointConstraints::New(pointDisplacement.mesh()).constrainDisplacement(pointDisplacement);
                    // tmp<pointField> ttransformedPts(new pointField(mesh.points()));
                    // pointField& transformedPts = ttransformedPts.ref();
          
                    // UIndirectList<point>(transformedPts, zm.pointIDs()) =
                    // pointField(newPoints.ref(), zm.pointIDs());
  
                    mesh.movePoints(newPoints);
                    //mesh.movePoints(mesh.points());

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

                // Solve the Momentum equation
                MRF.correctBoundaryVelocity(U);

                tmp<fvVectorMatrix> tUEqn
                (
                    fvm::ddt(U) + fvm::div(phi, U)
                    + MRF.DDt(U)
                    + turbulence->divDevReff(U)
                    ==
                    fvOptions(U)
                );
                fvVectorMatrix& UEqn = tUEqn.ref();

                UEqn.relax();
                fvOptions.constrain(UEqn);
                List<Eigen::MatrixXd> RedLinSysU;
                if (pimple.momentumPredictor())
                {
                    //solve(UEqn == -fvc::grad(p));
                    RedLinSysU = Umodes.project(UEqn, NmodesUproj);
                    volVectorField gradpfull = -fvc::grad(p);
                    Eigen::MatrixXd projGrad = Umodes.project(gradpfull, NmodesUproj);
                    RedLinSysU[1] = RedLinSysU[1] + projGrad;
                    a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual);
                    Umodes.reconstruct(U, a, "U");
                    fvOptions.correct(U);
                }

                // --- Pressure corrector loop
                while (pimple.correct())
                {

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
                    List<Eigen::MatrixXd> RedLinSysP;
                    bold = b;

                    // Non-orthogonal pressure corrector loop
                    while (pimple.correctNonOrthogonal())
                    {
                        fvScalarMatrix pEqn
                        (
                            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA) //p
                        );

                        pEqn.setReference(pRefCell, pRefValue);
                        //pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter()))); //p
                        RedLinSysP = Pmodes.project(pEqn, NmodesPproj);
                        b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
                        Pmodes.reconstruct(p, b, "p");

                        if (pimple.finalNonOrthogonalIter())
                        {
                            phi = phiHbyA - pEqn.flux();
                        }
                    }


                    // Explicitly relax pressure for momentum corrector
                    p.relax();
                    // b = bold + mesh.fieldRelaxationFactor("p")*(b-bold);
                    // Pmodes.reconstruct(p, b, "p");
                    U = HbyA - rAtU * fvc::grad(p); //p
                    U.correctBoundaryConditions();
                    fvOptions.correct(U);
                    // Correct Uf if the mesh is moving
                    fvc::correctUf(Uf, U, phi);

                    // Make the fluxes relative to the mesh motion
                    fvc::makeRelative(phi, U);

                }// end of the pimple.correct()
                Umodes.reconstruct(U, a, "U"); 
                Pmodes.reconstruct(p, b, "p");

            }// end of the pimple.loop()
            if(checkWrite(runTime))
            {
                romforcex.append(romforces.forceEff().x());
                romforcey.append(romforces.forceEff().y()); 

                centerofmassx.append(problem->sDRBMS().motion().centreOfMass().x());
                centerofmassy.append(problem->sDRBMS().motion().centreOfMass().y());
                centerofmassz.append(problem->sDRBMS().motion().centreOfMass().z());
                // To append the linear velocities
                velx.append(problem->sDRBMS().motion().v().x());
                vely.append(problem->sDRBMS().motion().v().y());
                velz.append(problem->sDRBMS().motion().v().z());

                ITHACAstream::exportSolution(U, name(counter), folder);
                ITHACAstream::exportSolution(p, name(counter), folder);
                ITHACAstream::exportSolution(pointDisplacement, name(counter), folder);
                ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
                std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
                UredFields.append(U.clone());
		        PredFields.append(p.clone());
		        Dfield.append(pointDisplacement.clone());
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

    void PodIpointDispl(Eigen::MatrixXd  coeffL2, Eigen::MatrixXd muu,  label NPdModes)
    {
        if (NPdModes == 0)
        {
            NPdModes = problem->Dmodes.size();
        }
        
        //coeffL2 = ITHACAutilities::getCoeffs(Dfield, Dmodes, NPdModes, false);

        //std::cout << "////////////////////////////////// :" << coeffL2 << std::endl;
        problem->samples.resize(NPdModes);
        problem->rbfSplines.resize(NPdModes);
        Eigen::MatrixXd weights;
      

        for (label i = 0; i < NPdModes; i++) // i is the nnumber of th mode
        {
            word weightName = "wRBF_M" + name(i + 1);

            if (ITHACAutilities::check_file("./ITHACAoutput/weights/" + weightName))
            {
                problem->samples[i] = new SPLINTER::DataTable(1, 1);
      //std::cout << "////////////////////////////////// :" << coeffL2.cols() << std::endl;
                for (label j = 0; j < coeffL2.cols();
                        j++) // j is the number of the nut snapshot
                {
                    problem->samples[i]->addSample(muu.row(j), coeffL2(i, j));
                }

                ITHACAstream::ReadDenseMatrix(weights, "./ITHACAoutput/weights/", weightName);
                problem->rbfSplines[i] = new SPLINTER::RBFSpline(*problem->samples[i],
                                                        SPLINTER::RadialBasisFunctionType::GAUSSIAN, weights);
                std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
            } 
           
            else
            {
               
                problem->samples[i] = new SPLINTER::DataTable(1, 1);

                for (label j = 0; j < coeffL2.cols();j++) // j is the number of the nut snapshot
                {
                    //std::cout << "///////////////////// :" << muu.row(j) << "///////////////////// :" << coeffL2(i, j) << std::endl;
                    problem->samples[i]->addSample(muu.row(j), coeffL2(i, j));
                }
                
                problem->rbfSplines[i] = new SPLINTER::RBFSpline(*problem->samples[i],
                                                        SPLINTER::RadialBasisFunctionType::MULTIQUADRIC);
                ITHACAstream::SaveDenseMatrix(problem->rbfSplines[i]->weights,
                                              "./ITHACAoutput/weights/", weightName);
            
                std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
            }
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

    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 50);
    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 50);
    int NmodesDout = para->ITHACAdict->lookupOrDefault<int>("NmodesDout", 10);

    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 25);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 25);
    int NmodesDproj = para->ITHACAdict->lookupOrDefault<int>("NmodesDproj", 5);



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
    example.finalTime = 10;
    example.timeStep = 0.01; //0.01;
    example.writeEvery = 0.1;

    // //Perform the offline solve
    example.offlineSolve();

    // // Read the lift functions
    //ITHACAstream::read_fields(example.liftfield, example.U, "./lift/");
    // //Search the lift function
    // example.liftSolve3();
    // // // Normalize the lifting function
    // ITHACAutilities::normalizeFields(example.liftfield);
    // // //Create homogeneous basis functions for velocity
    // online.computeLift(example.Ufield, example.liftfield, online.Uomfield);
    // // Perform a POD decomposition for velocity and pressure
    // ITHACAPOD::getModes(online.Uomfield, online.Umodes, online._U().name(),
    //                         example.podex, 0, 0, NmodesUout);

    //Perform POD on velocity pressure store the first 20 modes

    ITHACAPOD::getModes(example.Ufield, example.Umodes, online._U().name(),
                        example.podex, 0, 0, NmodesUout);


    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0,NmodesPout);
  
   ITHACAPOD::getModes(example.Dfield, example.Dmodes, online._pd().name(),
                        example.podex, 0, 0, NmodesDout);
   Eigen::MatrixXd FomCentersOfMassy = Foam2Eigen::field2Eigen(example.centerofmassy);
   //online.coeffL2 = ITHACAutilities::getCoeffs(example.Dfield, online.Dmodes, NmodesDproj, false);
   example.coeffL2 = ITHACAutilities::getCoeffs(example.Dfield, example.Dmodes, NmodesDproj, false);
   ITHACAstream::exportMatrix(online.coeffL2 , "coeffpd", "python","./ITHACAoutput/Matrices/");
   //online.PodIpointDispl(online.coeffL2, FomCentersOfMassy, NmodesDproj);
    Eigen::VectorXd FomVely = Foam2Eigen::field2Eigen(example.vely);
    //Eigen::VectorXd FomCentersOfMassy = Foam2Eigen::field2Eigen(example.centerofmassy);
    //online.CylDispl = FomCentersOfMassy;
    example.CylDispl = FomCentersOfMassy;
    //online.coeffL2 = 
    //reduced.reconstructLiftAndDrag(CoeffU, Coeffp, "./ITHACAoutput/forces");
    /// exporting matrices
    ITHACAstream::exportMatrix(FomVely, "FomVely", "python","./ITHACAoutput/VelOfDispl/");
    ITHACAstream::exportMatrix(FomCentersOfMassy, "FomCentersOfMassy", "python","./ITHACAoutput/CenterOfMass/");
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


    //Perform the online solutions
    //for (label k = 0; k < (example.mu).size(); k++)
    //{
    //scalar mu_now = example.mu(0, k);
    scalar mu_now = example.mu(0, 0);
    //example.change_viscosity(mu_now);
    //reduced.OnlineVelocity(vel);
    
    //example.restart();
    // ITHACAstream::writePoints(example.meshPtr().points(), "prova", "/polyMesh/");
    reduced.solveOnline_Pimple(mu_now, NmodesUproj, NmodesPproj, NmodesDproj);
    //}

    //scalar mu_now = example.mu(0, 0);
    //example.change_viscosity(mu_now);
    //reduced.OnlineVelocity(vel);
    //reduced.solveOnline_Pimple(mu_now, NmodesUproj, NmodesPproj);
    std::cout << "======================= ONLINE PHASE COMPLETED ================================" << "\n";
    /// Convert Foam List to Eigen Matrices
    Eigen::VectorXd RedVely = Foam2Eigen::field2Eigen(reduced.vely);
    Eigen::VectorXd RedCentersOfMassy = Foam2Eigen::field2Eigen(reduced.centerofmassy);
    
    ITHACAstream::exportMatrix(RedVely, "RedVely", "python","./ITHACAoutput/VelOfDispl/");
    ITHACAstream::exportMatrix(RedCentersOfMassy, "RedCentersOfMassy", "python","./ITHACAoutput/CenterOfMass/");
    
    std::cout << "======================= errorL2Rel ================================" << "\n";
    Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Rel(example.Ufield, reduced.UredFields);
    std::cout << "======================= errL2U completed================================" << "\n";
    Eigen::MatrixXd errL2P = ITHACAutilities::errorL2Rel(example.Pfield, reduced.PredFields);
    std::cout << "======================= errL2P completed================================" << "\n";

    Eigen::MatrixXd relErrorPdispl(example.Dfield.size(), 1);
    dimensionedVector pd("pd",dimLength,  vector(0, 1, 0));
    //dimensionedVector pd("pd",example._pd().dimensions(),  Zero);

    // for (label k = 0; k < example.Dfield.size(); k++)
    // {
    //   pointVectorField errorPdispl = (example.Dfield[k] - reduced.Dfield[k]).ref();
    //   pointVectorField devU = (example.Dfield[k] - pd).ref();
    //   relErrorPdispl(k, 0) = ITHACAutilities::frobNorm(errorPdispl)/
    //                        ITHACAutilities::frobNorm(devU);
    // }
    // ITHACAstream::exportMatrix(relErrorPdispl, "ErrorPdispl", "python","./ITHACAoutput/ErrorPdispl/");

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
