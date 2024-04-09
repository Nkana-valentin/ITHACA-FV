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
    Example of turbulent unsteady NS Reduction Problem solved by the use of the PIMPLE algorithm
SourceFiles
    PimpleRANSNNsRbf.C
\*---------------------------------------------------------------------------*/

// #include <torch/script.h>
// #include "torch2Eigen.H"
// #include "ITHACAstream.H"
// #include "ITHACAPOD.H"
// #include "fsiBasic.H"
//#include "ReducedBasicFsi.H"
#include "UnsteadyNSPimpleNN.H"

#include "pointConstraints.H"
#include "forces.H"
#include "IOmanip.H"
//#include "RBFMotionSolver.H"
#include <Eigen/SparseLU>
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
//#include <primitivePatchInterpolation.H>
#include <pointPatchField.H>
//#include <volPointInterpolation.H>
#include "quaternion.H"
#include <points0MotionSolver.H>


class reducedPimpleNNs
{
        public:
            explicit reducedPimpleNNs(UnsteadyNSPimpleNN& FoamPb): problem(&FoamPb)
                
            {

                   
                    for (int i = 0; i < problem->Umodes.size(); i++)
                    {
                        Umodes.append((problem->Umodes.toPtrList()[i]).clone());
                    }

                    for (int i = 0; i < problem->Pmodes.size(); i++)
                    {
                        Pmodes.append((problem->Pmodes.toPtrList()[i]).clone());
                    }

                    for (int i = 0; i < problem->nutModes.size(); i++)
                    {
                        nutModes.append((problem->nutModes.toPtrList()[i]).clone());
                    }

                  
                    for (int i = 0; i < problem->Dmodes.size(); i++)
                    {
                        Dmodes.append((problem->Dmodes.toPtrList()[i]).clone());
                    }
            }

            UnsteadyNSPimpleNN* problem;
            volScalarModes Pmodes;
            volVectorModes Umodes;
            pointVectorModes Dmodes;
            volScalarModes nutModes, NuTbar;

            
            ///////// time control variables
            scalar startTime = 0.0;
            scalar finalTime = 0.0;
            scalar timeStep = 0.0;
            scalar writeEvery = timeStep;
            scalar nextWrite = 0.0;
            /// List scalar for access the centerofmass
            List<scalar> centerofmassx, centerofmassy;
            //List<scalar> centerofmassz;
            List<scalar> Rx, Ry;
             /// List to save lift and drag forces
            List<scalar> romforcey, romforcex;
            // for saving pod interpolation  coefficients using rbf
            List<scalar> pdcoeffrbf;
            /// List of POD coefficients
            List<Eigen::MatrixXd> CoeffU, CoeffP, CoeffsNut; 
            PtrList<volScalarField> Pfield, nutFields;;
            PtrList<volVectorField> Ufield;
           
            PtrList<surfaceScalarField> phiFields;
            PtrList<pointVectorField> DfieldRbf;


          label counter = 1;
        void SolveOnlinePimple(scalar mu_now, 
                            int NmodesUproj, 
                            int NmodesPproj,
                            int NmodesNut,
                            int NmodesDproj, 
                            fileName folder = "./ITHACAoutput/Online/")
        {

            Eigen::MatrixXd a = Eigen::VectorXd::Zero(NmodesUproj);
            Eigen::MatrixXd b = Eigen::VectorXd::Zero(NmodesPproj);
            Eigen::MatrixXd nutCoeff = Eigen::VectorXd::Zero(NmodesNut);
            Time& runTime = problem->_runTime();
            dynamicFvMesh& mesh = problem->meshPtr();
            fv::options& fvOptions = problem->_fvOptions();
            pimpleControl& pimple = problem->_pimple();
            volVectorField& U = problem->_U();
            volScalarField& p = problem->_p();
            surfaceScalarField& phi = problem->_phi(); 
            pointVectorField& pointDisplacement = problem->_pointDisplacement();
            IOMRFZoneList& MRF = problem->_MRF();
            singlePhaseTransportModel& laminarTransport = problem->_laminarTransport();
            autoPtr<incompressible::turbulenceModel> turbulence(incompressible::turbulenceModel::New(U, 
                phi, laminarTransport));
            volScalarField& nut = problem->_nut();
            volScalarField nueff = turbulence->nuEff();
            instantList Times = runTime.times();
            runTime.setEndTime(finalTime);

// To solve the RBFI system to obtain the weights
//PodIpointDispl(problem->coeffL2, problem->CylDispl, problem->CylRot, NmodesDproj);

            runTime.setTime(Times[1], 1);
            runTime.setDeltaT(timeStep);
            nextWrite = startTime;
            label pRefCell = 0;
            scalar pRefValue = 0.0;

            bool  correctPhi = problem->correctPhi;
            bool  checkMeshCourantNo = problem->checkMeshCourantNo;
            bool  moveMeshOuterCorrectors = problem->moveMeshOuterCorrectors;
            scalar  cumulativeContErr = problem->cumulativeContErr;
#include "createUfIfPresent.H"
            turbulence->validate();
            dictionary dictCoeffs(problem->dyndict->findDict("sixDoFRigidBodyMotionCoeffs")); 
            sixDoFRigidBodyMotion sDRBM(dictCoeffs, dictCoeffs, runTime);
            vector rotationAngle(Zero);

            Foam::functionObjects::forces romforces("romforces", mesh,  dictCoeffs);
            Eigen::VectorXd  pdCoeffNew;
            pdCoeffNew.resize(NmodesDproj);
            label curTimeIndex_ = -1;

            // PIMPLE algorithm starts here
            Info<< "\nStarting time loop\n" << endl;
            while (runTime.run())
            {
//#include "readDyMControls.H"
#include "CourantNo.H"

                runTime.setEndTime(finalTime);
                runTime++;
                //p.storePrevIter();
                Info << "Time = " << runTime.timeName() << nl << endl; 
                /// nuEff calculation.
                nueff = nut + turbulence->nu();

                while (pimple.loop())
                {
                    if (pimple.firstIter() || moveMeshOuterCorrectors)
                    {
            
//#include "test.H"
                            mesh.controlledUpdate();
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

                        fvVectorMatrix UEqn(
                        fvm::ddt(U) 
                        + fvm::div(phi, U) 
                        //+ MRF.DDt(U)
                        //+ turbulence->divDevReff(U) 
                        - fvm::laplacian(nueff, U)
                        - fvc::div(nueff * dev2(T(fvc::grad(U))))
                        ==  
                            fvOptions(U)
                    );
                    //fvVectorMatrix& UEqn = tUEqn.ref();
                    UEqn.relax();
                    fvOptions.constrain(UEqn);

                    List<Eigen::MatrixXd> RedLinSysU;
                    if (pimple.momentumPredictor())
                    {

                        RedLinSysU = Umodes.project(UEqn, NmodesUproj, "G");///???
                        volVectorField gradpfull = -fvc::grad(p);
                        Eigen::MatrixXd projGrad = Umodes.project(gradpfull, NmodesUproj);
                        RedLinSysU[1] = RedLinSysU[1] + projGrad;
                        a = RedLinSysU[0].householderQr().solve(RedLinSysU[1]);
                        //a = RedLinSysU[0].colPivHouseholderQr().solve(RedLinSysU[1]);
                        //volVectorField U("U", Ubar[0]);
                        Umodes.reconstruct(U, a, "U");
                        //solve(UEqn == -fvc::grad(p));
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

                        // Update the pressure BCs to ensure flux consistency
                        constrainPressure(p, U, phiHbyA, rAtU(), MRF); //p
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
                            RedLinSysP = Pmodes.project(pEqn, NmodesPproj, "G" );
                            /// Solve for the reduced coefficient for pressure
                            b = RedLinSysP[0].householderQr().solve(RedLinSysP[1]); 
                            //b = RedLinSysP[0].colPivHouseholderQr().solve(RedLinSysP[1]);
                            Pmodes.reconstruct(p, b, "p");

                           if (pimple.finalNonOrthogonalIter())
                            {
                                phi = phiHbyA - pEqn.flux();
                            }
                        }
                        // Explicitly relax pressure for momentum corrector
                        p.relax();
                        U = HbyA - rAtU * fvc::grad(p); //p
                        U.correctBoundaryConditions();
                        fvOptions.correct(U);
                        // Correct Uf if the mesh is moving
                        fvc::correctUf(Uf, U, phi);
                        // Make the fluxes relative to the mesh motion
                        fvc::makeRelative(phi, U);

                    }// end of the pimple.correct()

                }// end of the pimple.loop()
                /// Eval the Lstm model                     
                /// Nut reconstruction using NNs
                nutCoeff = problem->evalNet(a);
                // // nutCoeff = problem->evalLstm(nutCoeff);
                nutModes.reconstruct(nut, nutCoeff, "nut");
                // //nut = nut + NuTbar[0];
                // nutCoeff = ITHACAutilities::getCoeffs(nut, nutModes, NmodesNut, true);
                //nut =  turbulence->nut();
                if(checkWrite(runTime))
                {
                    
                    Rx.append( sDRBM.centreOfRotation().y() );
                    Ry.append( rotationAngle.z() );
                    centerofmassy.append(sDRBM.centreOfMass().y() );
                    
                    // pdcoeffrbf.append(pdCoeff(0, 0));
                    CoeffP.append(b);
                    CoeffU.append(a);
                    CoeffsNut.append(nutCoeff);
                    romforcey.append(romforces.forceEff().y());
                    romforcex.append(romforces.forceEff().x());

                    ITHACAstream::exportSolution(U, name(counter), folder);
                    ITHACAstream::exportSolution(p, name(counter), folder);
                    //ITHACAstream::exportSolution(phi, name(counter), folder);
                    ITHACAstream::exportSolution(pointDisplacement, name(counter), folder);
                    ITHACAstream::exportSolution(nut, name(counter), folder);
                    ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
                    std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
                    Ufield.append(U.clone());
                    Pfield.append(p.clone());
                    nutFields.append(nut.clone() );
                    //phiFields.append(phi.clone());
                    DfieldRbf.append(pointDisplacement.clone() );
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

       
        void PodIpointDispl(Eigen::MatrixXd coeffL2,Eigen::MatrixXd mu1, Eigen::MatrixXd mu2, label NPdModes)
        {


            if (NPdModes == 0)
            {
                NPdModes = Dmodes.size();
            }

            problem->samples.resize(NPdModes);
            problem->rbfSplines.resize(NPdModes);
            Eigen::MatrixXd weights;
            std::vector<double> x(2);
          
            for (label i = 0; i < NPdModes; i++) // i is the nnumber of th mode
            {
                word weightName = "wRBF_M" + name(i + 1);

                if (ITHACAutilities::check_file("./ITHACAoutput/weights/" + weightName))
                {
                    problem->samples[i] = new SPLINTER::DataTable(true, true);
                    for (label j = 0; j < coeffL2.cols(); j++) // j is the number of the nut snapshot
                    {

                        x[0] = mu1(j);
                        x[1] = mu2(j);
                        problem->samples[i]->addSample(x, coeffL2(i, j));
                        //problem->samples[i]->addSample(mu1.row(j), coeffL2(i, j));
                    }

                    ITHACAstream::ReadDenseMatrix(weights, "./ITHACAoutput/weights/", weightName);
                    problem->rbfSplines[i] = new SPLINTER::RBFSpline(*(problem->samples)[i], SPLINTER::RadialBasisFunctionType::THIN_PLATE_SPLINE, weights);
                    std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
                } 
               
                else
                {
                    problem->samples[i] = new SPLINTER::DataTable(true, true);

                    for (label j = 0; j < coeffL2.cols();j++) // j is the number of the nut snapshot
                    {
                        x[0] = mu1(j);
                        x[1] = mu2(j);
                        problem->samples[i]->addSample(x, coeffL2(i, j));
                        //problem->samples[i]->addSample(mu1.row(j), coeffL2(i, j));
                    }
                    problem->rbfSplines[i] = new SPLINTER::RBFSpline(*(problem->samples)[i], SPLINTER::RadialBasisFunctionType::THIN_PLATE_SPLINE);
                    ITHACAstream::SaveDenseMatrix(problem->rbfSplines[i]->weights, "./ITHACAoutput/weights/", weightName);
                    std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
                }
            }
        }

};            


class unsteadypimplernns : public UnsteadyNSPimpleNN
{
    public:
        /// Constructor
        explicit unsteadypimplernns(int argc, char* argv[])
            :
            UnsteadyNSPimpleNN(argc, argv),
            U(_U()),p(_p()),
            nut(_nut()), 
            pD(_pointDisplacement())
           {
            // curX = _mesh().points();
            // point0 = _mesh().points();
            //std::cout << "################ line 789 ##################" << std::endl;
          }
            // Relevant Fields
         volVectorField& U;
         volScalarField& p;
         volScalarField& nut;
         pointVectorField& pD;

        /// Perform an Offline solve
        void offlineSolve(word folder="./ITHACAoutput/Offline/")
        {
                List<scalar> mu_now(1);

                if (offline)
                {
                    ITHACAstream::read_fields(Ufield, U, folder );
                    ITHACAstream::read_fields(Pfield, p, folder );
                    ITHACAstream::read_fields(nutFields, nut, folder );
                    ITHACAstream::read_fields(Dfield, pD, folder);
                    //ITHACAstream::read_fields(ListPointsField, DataPoints(), folder);
                     // mu_samples =
                     // ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
                }
                else
                {
                    mu_now[0] = mu(0, 0); //mu.cols()=50
                    truthSolve(mu_now, folder);
                    //restart();
                }    

        }

};


int main(int argc, char* argv[])
{
    // Construct the tutorial object
    unsteadypimplernns example(argc, argv);
    // Info << example._U() << endl;
    // Info << example._p() << endl;
    //return 0; 
    unsteadypimplernns online(argc,argv);
    std::clock_t startOff;
    double durationOff;
    // Read some parameters from file
    ITHACAparameters* para = ITHACAparameters::getInstance(example.meshPtr(),
                                                            example._runTime());
    int NmodesUout = para->ITHACAdict->lookupOrDefault<int>("NmodesUout", 50);
    int NmodesUproj = para->ITHACAdict->lookupOrDefault<int>("NmodesUproj", 10);

    int NmodesPout = para->ITHACAdict->lookupOrDefault<int>("NmodesPout", 50);
    int NmodesPproj = para->ITHACAdict->lookupOrDefault<int>("NmodesPproj", 10);

    int NmodesNutOut = para->ITHACAdict->lookupOrDefault<int>("NmodesNutOut", 7);
    int NmodesNutProj = para->ITHACAdict->lookupOrDefault<int>("NmodesNutProj", 7);
    // std::cout << NmodesPproj << std::endl;
    // std::cout << NmodesUproj << std::endl;
    // exit(0);
    int NmodesDout = para->ITHACAdict->lookupOrDefault<int>("NmodesDout", 10);
    int NmodesDproj = para->ITHACAdict->lookupOrDefault<int>("NmodesDproj", 3);

    // Read the par file where the parameters are stored
    word filename("./par");
    example.mu = ITHACAstream::readMatrix(filename);

    //Set the inlet boundaries patch 0 directions x and y
    example.inletIndex.resize(1, 2);
    example.inletIndex(0, 0) = 0;
    example.inletIndex(0, 1) = 0;
    // Time parameters: We can use Ioodictionnary to access time parameters
    example.startTime  = 0;
    example.finalTime  = 1;//0.5; //0.5; //0.05; //0.5;  
    example.timeStep   = 5e-06; //2e-06
    example.writeEvery = 5e-04; //0.00001; //1e-02; //5e-04, 5e-03
    // //Perform the offline solve
    startOff= std::clock();
    example.offlineSolve();
    //exit(0);
    durationOff = (std::clock() - startOff);
    std::cout << "The Offline phase  duration  is  equal  to " << durationOff << std::endl;

    if(!ITHACAutilities::check_folder("./ITHACAoutput/DataFromFoam"))
    {

        mkDir("ITHACAoutput/DataFromFoam");
        Eigen::VectorXd fomforcex = Foam2Eigen::field2Eigen(example.fomforcex);
        //ITHACAstream::exportMatrix(fomforcex, "fomforcex", "python","./ITHACAoutput/DataFromFoam/");
        cnpy::save(fomforcex, "./ITHACAoutput/DataFromFoam/fomforcex.npy");

        Eigen::VectorXd fomforcey = Foam2Eigen::field2Eigen(example.fomforcey);
        //ITHACAstream::exportMatrix(fomforcey, "fomforcey", "python","./ITHACAoutput/DataFromFoam/");
        cnpy::save(fomforcey, "./ITHACAoutput/DataFromFoam/fomforcey.npy");

        Eigen::MatrixXd Ry = Foam2Eigen::field2Eigen(example.Ry);
        //ITHACAstream::exportMatrix(Ry, "Ry", "python","./ITHACAoutput/DataFromFoam/");
        cnpy::save(Ry, "./ITHACAoutput/DataFromFoam/Ry.npy");

        Eigen::MatrixXd Rx = Foam2Eigen::field2Eigen(example.Rx);
        //ITHACAstream::exportMatrix(Rx, "Rx", "python","./ITHACAoutput/DataFromFoam/");
        cnpy::save(Rx, "./ITHACAoutput/DataFromFoam/Rx.npy");

        Eigen::MatrixXd CentreOfMassY = Foam2Eigen::field2Eigen(example.centerofmassy);
        //ITHACAstream::exportMatrix(CentreOfMassY, "CentreOfMassY", "python", "./ITHACAoutput/DataFromFoam");
        cnpy::save(CentreOfMassY, "./ITHACAoutput/DataFromFoam/CentreOfMassY.npy");
    }

    if(std::ifstream("./ITHACAoutput/DataFromFoam/Ry.npy"))
    {
        Info << "################ Reading pitch ##############" << endl;
        cnpy::load(online.CylRot, "./ITHACAoutput/DataFromFoam/Ry.npy");
    }

    if(std::ifstream("./ITHACAoutput/DataFromFoam/Rx.npy"))
    {
        Info << "################ Reading plunge ##############" << endl;
        cnpy::load(online.CylDispl, "./ITHACAoutput/DataFromFoam/Rx.npy");
    }
   
    int const p = 1001;
    int const q = 1;

    int i = 0;
    int j = 0;

    int s = 0;
    int e = 1000;

    Eigen::MatrixXd plunge = online.CylDispl.block<p,q>(i,j); //online.CylDispl.topRows<500>();
    Eigen::MatrixXd pitch  = online.CylRot.block<p,q>(i,j);   // online.CylRot.topRows<500>();

    /// Update the plunge and pitch matrices 
    online.CylDispl = plunge;
    online.CylRot   = pitch;

    // Info << "plunge.size() = \t" <<  plunge.size() << endl;
    // Info << "pitch.size() = \t" <<   pitch.size() << endl;

    //std::cout << "online.CylDispl = \n" << online.CylDispl  << std::endl;
    // std::cout << "plunge = \t" <<   online.CylDispl << std::endl;

    // //std::cout << online.CylDispl.cols() << std::endl;

    //exit(0);
    // volScalarField NuTbar = ITHACAutilities::computeAverage(example.nutFields);
    // PtrList<volScalarField> NuTcentre;
    // for (label i = 0; i < example.nutFields.size(); i++){
    //         volScalarField nut("nut", example.nutFields[i] - NuTbar );
    //         NuTcentre.append(nut.clone() );
    // }

    PtrList<volScalarField> Pfield, nutFields;
    PtrList<volVectorField> Ufield;
    PtrList<pointVectorField> Dfield;

    for (label i = s; i < e+1; i++)
    {
        Ufield.append(example.Ufield[i].clone() );
        Pfield.append(example.Pfield[i].clone() );
        nutFields.append(example.nutFields[i].clone() );
        Dfield.append(example.Dfield[i].clone() );
    }

    // Info << "Ufield.size() = \t" <<  Ufield.size() << endl;
    // Info << "Pfield.size() = \t" <<  Pfield.size() << endl;
    // Info << "nutFields.size() = \t" <<  nutFields.size() << endl;
    // Info << "Dfield.size() = \t" <<  Dfield.size() << endl;


    // exit(0);

    if(!ITHACAutilities::check_folder("./ITHACAoutput/POD/1"))
    {

        ITHACAPOD::getModes(Ufield, online.Umodes, example._U().name(),
                                    example.podex, 0, 0, NmodesUout);
        ITHACAPOD::getModes(Pfield, online.Pmodes, example._p().name(),
                                    example.podex, 0, 0,NmodesPout);
        ITHACAPOD::getModes(nutFields, online.nutModes, example._nut().name(),
                                    example.podex, 0, 0, NmodesNutOut); 
       // ITHACAPOD::getModes(Dfield, online.Dmodes, online._pointDisplacement().name(),
       //                 example.podex, 0, 0, NmodesDout);
    }

    /// Create the NNs for nut
    online.getTurbNN();
    /// Get temporal coefficients for pointDisplacement fields
    //online.coeffL2 = ITHACAutilities::getCoeffs(Dfield, online.Dmodes, NmodesDproj, false);

    // exit(0);

    // if(!ITHACAutilities::check_folder("./ITHACAoutput/POD/1"))
    // {

    //     ITHACAPOD::getModes(example.Ufield, online.Umodes, example._U().name(),
    //                                 example.podex, 0, 0, NmodesUout);
    //     ITHACAPOD::getModes(example.Pfield, online.Pmodes, example._p().name(),
    //                                 example.podex, 0, 0,NmodesPout);
    //     ITHACAPOD::getModes(example.nutFields, online.nutModes, example._nut().name(),
    //                                 example.podex, 0, 0, NmodesNutOut); 
    //    ITHACAPOD::getModes(example.Dfield, online.Dmodes, online._pointDisplacement().name(),
    //                    example.podex, 0, 0, NmodesDout);
    // }

    // /// Create the NNs for nut
    // online.getTurbNN();
    // /// Get temporal coefficients for pointDisplacement fields
    // online.coeffL2 = ITHACAutilities::getCoeffs(example.Dfield, online.Dmodes, NmodesDproj, false);
    /// Load the network by training the model.
    online.loadNet("ITHACAoutput/NN/Net_" + name(example.NUmodes) + "_" + name(example.NNutModes) + ".pt");
    exit(0);
    //online.loadLstmNet("ITHACAoutput/NN/LstmNet_" + name(example.NNutModes) + ".pt");

    /// ############### contruct the reduced the class object ###################
    reducedPimpleNNs reduced(online);
    //reduced.NuTbar.append(NuTbar.clone());

    reduced.startTime  =  example.startTime;
    reduced.finalTime  =  example.finalTime; //+1
    reduced.timeStep   =  example.timeStep;
    reduced.writeEvery = example.writeEvery;
    scalar mu_now = example.mu(0, 0);

    std::clock_t startOn;
    double durationOn;
    startOn = std::clock();


    reduced.SolveOnlinePimple(mu_now, NmodesUproj, NmodesPproj, example.NNutModes,  NmodesDproj);
    durationOn = (std::clock() - startOn);
    std::cout << "The Online  phase  duration  is  equal  to " << durationOn << std::endl;

    if (!ITHACAutilities::check_folder("./ITHACAoutput/DataFromRom"))
    {
            mkDir("ITHACAoutput/DataFromRom");

    }
    Eigen::VectorXd romforcex = Foam2Eigen::field2Eigen(reduced.romforcex);
    cnpy::save(romforcex, "./ITHACAoutput/DataFromRom/romforcex_" + name(example.NUmodes) + "_" + name(example.NPmodes) + ".npy");

    Eigen::VectorXd romforcey = Foam2Eigen::field2Eigen(reduced.romforcey);
    cnpy::save(romforcey, "./ITHACAoutput/DataFromRom/romforcey_" + name(example.NUmodes) + "_" + name(example.NPmodes) + ".npy");

    Eigen::VectorXd Rx = Foam2Eigen::field2Eigen(reduced.Rx);
    cnpy::save(Rx, "./ITHACAoutput/DataFromRom/Rx_" + name(example.NUmodes) + "_" + name(example.NPmodes) + ".npy");

    Eigen::VectorXd Ry = Foam2Eigen::field2Eigen(reduced.Ry);
    cnpy::save(Ry, "./ITHACAoutput/DataFromRom/Ry_" + name(example.NUmodes) + "_" + name(example.NPmodes) + ".npy");

    Eigen::MatrixXd RedCentersOfMassy = Foam2Eigen::field2Eigen(reduced.centerofmassy);
    cnpy::save(RedCentersOfMassy, "./ITHACAoutput/DataFromRom/RedCentersOfMassy_" + name(example.NUmodes) + "_" + name(example.NPmodes) + ".npy");
    /// Export matrix for reduced coeffs
    ITHACAstream::exportMatrix(reduced.CoeffsNut,"CoeffsNut", "python", "./ITHACAoutput/DataFromRom/");
    //cnpy::save(reduced.CoeffsNut, "./ITHACAoutput/DataFromRom/CoeffsNut.npy");
    ITHACAstream::exportMatrix(reduced.CoeffU,   "CoeffU",    "python", "./ITHACAoutput/DataFromRom/");
    //cnpy::save(reduced.CoeffU, "./ITHACAoutput/DataFromRom/CoeffU.npy");
    ITHACAstream::exportMatrix(reduced.CoeffP,   "CoeffP",    "python", "./ITHACAoutput/DataFromRom/");
    //cnpy::save(reduced.CoeffP, "./ITHACAoutput/DataFromRom/CoeffP.npy");

    // //exit(0);
    //std::cout << "======================= L2 Relative Errors ================================" << "\n";
    Eigen::MatrixXd errL2U = ITHACAutilities::errorL2Rel(example.Ufield, reduced.Ufield);
    cnpy::save(errL2U, "./ITHACAoutput/DataFromRom/errL2U_" + name(example.NUmodes) + "_" + name(example.NPmodes) + "_" + name(example.NNutModes) + ".npy");

    std::cout << "======================= errL2U completed================================" << "\n";
    Eigen::MatrixXd errL2P = ITHACAutilities::errorL2Rel(example.Pfield, reduced.Pfield);
    cnpy::save(errL2P, "./ITHACAoutput/DataFromRom/errL2P_" + name(example.NUmodes) + "_" + name(example.NPmodes) + "_" + name(example.NNutModes) + ".npy");

    std::cout << "======================= errL2P completed================================" << "\n";
    Eigen::MatrixXd errL2Nut = ITHACAutilities::errorL2Rel(example.nutFields, reduced.nutFields);
    cnpy::save(errL2Nut, "./ITHACAoutput/DataFromRom/errL2Nut_" + name(example.NUmodes) + "_" + name(example.NPmodes)+ "_" + name(example.NNutModes) + ".npy");

    exit(0);
}
