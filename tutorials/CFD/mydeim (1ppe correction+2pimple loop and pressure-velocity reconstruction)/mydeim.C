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
    Example of the reconstruction of a non-linear function using the DEIM
SourceFiles
    FSIDEIM.C
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "sixDoFRigidBodyMotionSolver.H"
#include "ITHACAutilities.H"
#include "Foam2Eigen.H"
#include "DEIM.H"
#include "ITHACAstream.H"
#include <chrono>
#include "Modes.H"
#include "reductionProblem.H"
#include "ReducedProblem.H"
#include "EigenFunctions.H"
#include <Eigen/SVD>
#include <Eigen/SparseLU>
#include <Eigen/Dense>
#include <cmath>


class vVf : public DEIM<volVectorField>
{
    public:
        using DEIM::DEIM;
        autoPtr<volVectorField> field;
        autoPtr<volVectorField> subfield;

        Eigen::VectorXd onlineCoeffs(volVectorField& f)
        {
            theta.resize(magicPoints().size()); //Eigen::VectorXd::Zero(NmodesUproj);
            Eigen::MatrixXd ftoeigen = Foam2Eigen::field2Eigen(f);
            // for (int i = 0; i < magicPoints().size(); i++)
            // {
            //     Info << "===============" << f[localMagicPoints[i]] << "======================" <<endl;
            //     //Eigen::VectorXd ftoeigen = Foam2Eigen::field2Eigen(f[localMagicPoints[i]]);
            //     // std::cout << ftoeigen;
            // }
            theta = P.transpose() * ftoeigen;
            
            // Info << "==============================================================" <<endl;
            // std::cout << "theta rows "<< theta.rows() << std::endl;
            // std::cout << "theta cols "<< theta.cols() << std::endl;
            // Info << "==============================================================" <<endl;
            // std::cout << "ftoeigen rows "<< ftoeigen.rows() << std::endl;
            // std::cout << "ftoeigen cols "<< ftoeigen.cols() << std::endl;
            // Info << "==============================================================" <<endl;
            // std::cout <<"P rows "<<  P.rows() << std::endl;
            // std::cout <<"P cols "<<  P.cols() << std::endl;

            return theta;
        }

        void save_magic()
        {
            Eigen::VectorXd magic_indices;
            magic_indices.resize(magicPoints().size());
            for (label i = 0; i < magicPoints().size(); i++)
                magic_indices[i]=magicPoints()[i];

            ITHACAstream::exportMatrix(magic_indices, "magicPoint", "python",
                                           "./ITHACAoutput/magicPoints/");
        }


        
};
class DeimScalarField: public DEIM<volScalarField>
{
    public:
        using DEIM::DEIM;
        autoPtr<volScalarField> field;
        autoPtr<volScalarField> subfield;
        
        
};


class FSIDEIM : public reductionProblem
{
public:
        explicit FSIDEIM(int argc, char* argv[])
        {

            _args = autoPtr<argList>
                (
                    new argList(argc, argv)
                );

            if (!_args->checkRootCase())
            {
                Foam::FatalError.exit();
            }

            argList& args = _args();

#include "createTime.H"
              meshPtr = autoPtr<dynamicFvMesh> (dynamicFvMesh::New(args, runTime));
              dynamicFvMesh& mesh = meshPtr();
             _pimple = autoPtr<pimpleControl>
                   (
                       new pimpleControl
                       (
                           mesh
                       )
                   );
                ITHACAdict = new IOdictionary
                (
                    IOobject
                    (
                        "ITHACAdict",
                        runTime.system(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );
#include "createFields.H" 


            para = ITHACAparameters::getInstance(mesh, runTime); 
            point0 = mesh.points(); 
            NmodesU  = readInt(ITHACAdict->lookup("NmodesU"));
            NmodesS  = readInt(ITHACAdict->lookup("NmodesS")); 
            NmodesP  = readInt(ITHACAdict->lookup("NmodesP")); 

            std::cout << "####### test of the fsideim ctor ######" << std::endl;
        }
        //FSIDEIM() {};
        ~FSIDEIM() {

             std::cout << "####### fsideim dtor ######" << std::endl;
        };

         ///////// time control variables
        scalar startTime = 0.0;
        scalar finalTime = 0.0;
        scalar timeStep = 0.0;
        scalar writeEvery = timeStep;
        scalar nextWrite = 0.0;

        autoPtr<argList> _args;
        autoPtr<pimpleControl> _pimple;
        autoPtr<fv::options> _fvOptions;
        autoPtr<Time> _runTime;
        // Offline fields
        PtrList<volVectorField> Ufield;  
        PtrList<volVectorField> UlinearField; 
        PtrList<volVectorField> Tfield; 
        PtrList<volScalarField> Pfield;
        
        autoPtr<IOMRFZoneList> _MRF;
        autoPtr<volScalarField> _p;
        autoPtr<volVectorField> _U;
        autoPtr<volVectorField> _S;
        autoPtr<volVectorField> _T;
        autoPtr<pointVectorField> _pd;
        autoPtr<surfaceScalarField> _phi;
        autoPtr<surfaceScalarField> _phi0;
        autoPtr<Foam::dynamicFvMesh> meshPtr;
        bool  correctPhi;
        bool  checkMeshCourantNo;
        bool  moveMeshOuterCorrectors;
        ITHACAparameters* para;
        vectorField point0;
      
        autoPtr<incompressible::turbulenceModel> turbulence;
        autoPtr<singlePhaseTransportModel> _laminarTransport;
        autoPtr<sixDoFRigidBodyMotionSolver> sDRBMS;
        IOdictionary* dyndict;

        Eigen::MatrixXd ModesUEig;
        Eigen::MatrixXd ModesTEig;
        vVf* DEIMObject;
        vVf* Lt;
        DeimScalarField* DEIMPField;
        Eigen::MatrixXd ReducedVectors;

        volVectorModes Tmodes;
        volScalarModes Pmodes;
        // List of Modes for U
        volVectorModes Umodes;
        // List of Modes for Phi
        surfaceScalarModes Phimodes;

        int NmodesU;
        int NmodesS;
        int NmodesP;
        void offlineSolve(fileName folder = "./ITHACAoutput/Offline") 
        {

            Time& runTime = _runTime();
            surfaceScalarField& phi = _phi();
            dynamicFvMesh& mesh = meshPtr();
            fv::options& fvOptions = _fvOptions();
            pimpleControl& pimple = _pimple();
            volScalarField& p = _p();
            volVectorField& U = _U();
            volVectorField& T = _T();
            T = -fvc::div(_phi(), _U()); // -fvc::grad(_p());
           
            IOMRFZoneList& MRF = _MRF();
            singlePhaseTransportModel& laminarTransport = _laminarTransport();
            turbulence = autoPtr<incompressible::turbulenceModel>(incompressible::turbulenceModel::New(U, phi, laminarTransport));
            volScalarField nu = laminarTransport.nu();
            label pRefCell = 0;
            scalar pRefValue = 0.0;
            instantList Times = runTime.times();
            runTime.setEndTime(finalTime);
            runTime.setTime(Times[1], 1);
            runTime.setDeltaT(timeStep);
            nextWrite = startTime; 
            label counter = 1;
#include "createUfIfPresent.H"
            turbulence->validate();
            /// Exporting initial solution
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(T, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter ), folder);
            
            ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
            std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
            Ufield.append(U.clone());
            Tfield.append(T.clone());
            Pfield.append(p.clone());
            //Phifield.append(phi.clone());
            scalar  cumulativeContErr = 0.0;
            counter++;
            nextWrite += writeEvery;
  
            while(runTime.run()){

                runTime.setEndTime(finalTime);
                runTime++;
                //p.storePrevIter();
                Info<< "Time = " << runTime.timeName() << nl << endl;
       
                // --- Pressure-velocity PIMPLE corrector loop
                while (pimple.loop())
                {
                    #include "UEqn.H"

                    // --- Pressure corrector loop
                    while (pimple.correct())
                    {
                        #include "pEqn.H"
                    }

                    if (pimple.turbCorr())
                    {
                        laminarTransport.correct();
                        turbulence->correct();
                    }

                } // End pimple loop
                // Eval the non-linear term
                T = -fvc::div(phi, U); // - fvc::grad(p);        
                if(checkWrite(runTime))
                {
                    ITHACAstream::exportSolution(U, name(counter), folder);
                    ITHACAstream::exportSolution(T, name(counter), folder);
                    ITHACAstream::exportSolution(p, name(counter), folder);
                    ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
                    std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
                    Ufield.append(U.clone());
                    Pfield.append(p.clone());
                    Tfield.append(T.clone());
                    counter++;
                    nextWrite += writeEvery;
                }        
               
               
            }// end runTime loop

        }
  
        void PODDEIM(int NmodesS)
        {
            dynamicFvMesh& mesh = meshPtr();
            DEIMObject = new vVf(Tfield, NmodesS, "convectiveterm", _T().name());
            DEIMObject->subfield = autoPtr<volVectorField>(new volVectorField(DEIMObject->generateSubmesh(0, mesh, _T())));
           
        }
       

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

        void restart()
        {
                turbulence.clear();
                _fvOptions.clear();
                _laminarTransport.clear();
                _p.clear();
                _U.clear();
                _phi.clear();
                //_Uf.clear();
                _pd.clear();
                _T.clear();
                sDRBMS.clear();
                argList& args = _args();
                Time& runTime = _runTime();
                runTime.setTime(0, 1);
                //meshPtr().movePoints(point0);    
                //pointField& pointOld = const_cast<pointField&> (meshPtr().oldPoints());
                //pointOld = point0;
                // _pimple.clear();
                Foam::dynamicFvMesh& mesh = meshPtr();
 //std::cout << "=============================================================" << std::endl;
                _pimple = autoPtr<pimpleControl>
                          (
                                   new pimpleControl
                                   (
                                       mesh
                                   )
                           );

                
#include "createFields.H" 
        }
 

};

class ReducedFSIDEIM : public reductionProblem
{
        public:
                explicit ReducedFSIDEIM(FSIDEIM& FoamProblem) : problem(&FoamProblem)
                {

           std::cout << "============================       inside reduced deim=================================" << std::endl;
                    //for (int i = 0; i < problem->Umodes.size(); i++)
                    //{
                    //    Umodes.append((problem->Umodes.toPtrList()[i]).clone());
                    //}

                    ///for (int i = 0; i < problem->Pmodes.size(); i++)
                    //{
                    //    Pmodes.append((problem->Pmodes.toPtrList()[i]).clone());
                    //}
                    // for (int i = 0; i < problem->DEIMObject->modes.size(); i++)
                    // {
                    //     Tmodes.append((problem->Tmodes.toPtrList()[i]).clone());
                    // }
                     std::cout << "=============================================================" << std::endl;

                }
                // Pointer to the full problem
                FSIDEIM* problem;
                // Online fields
                PtrList<volVectorField> Ured;  
                PtrList<volScalarField> Pred;
                PtrList<volVectorField> Tred;
                volScalarModes Pmodes;
                volVectorModes Umodes;
                volVectorModes Tmodes;
                scalar startTime = 0.0;
                scalar finalTime = 0.0;
                scalar timeStep = 0.0;
                scalar writeEvery = timeStep;
                scalar nextWrite = 0.0;
                label counter = problem->counter;
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
                void onlineSolve(int NmodesU,  int NmodesP, int NmodesS, fileName folder = "./ITHACAoutput/Online") 
                {

                    Time& runTime = problem->_runTime();
                    surfaceScalarField& phi = problem->_phi();
                    dynamicFvMesh& mesh = problem->meshPtr();
                    fv::options& fvOptions = problem->_fvOptions();
                    pimpleControl& pimple = problem->_pimple();
                    volScalarField& p = problem->_p();
                    volVectorField& U = problem->_U();
                    volVectorField& T = problem->_T();
                    T = -fvc::div(problem->_phi(), problem->_U()); //- fvc::grad(problem->_p());
                    IOMRFZoneList& MRF = problem->_MRF();
                    singlePhaseTransportModel& laminarTransport = problem->_laminarTransport();
                    label pRefCell = 0;
                    scalar pRefValue = 0.0;
                    instantList Times = runTime.times();
                    runTime.setEndTime(finalTime);
                    runTime.setTime(Times[1], 1);
                    runTime.setDeltaT(timeStep);
                    nextWrite = startTime; 
                    label counter = 1;
                    Eigen::MatrixXd a = Eigen::VectorXd::Zero(NmodesU);
                    Eigen::MatrixXd b = Eigen::VectorXd::Zero(NmodesP);
                    Eigen::VectorXd uresidual = Eigen::VectorXd::Zero(NmodesU);
                    a = ITHACAutilities::getCoeffs(U, problem->Umodes, NmodesU, true);
                    Eigen::VectorXd presidual = Eigen::VectorXd::Zero(NmodesP);
                    b = ITHACAutilities::getCoeffs(p, problem->Pmodes, NmodesP, true);
                    Info << "###########################################################" << endl;
                    problem->Umodes.reconstruct(U, a, "U");
                    problem->Pmodes.reconstruct(p, b, "p");
                    
                    // Online evaluation of the non linear function
                    Eigen::VectorXd thetaon = Eigen::VectorXd::Zero(NmodesS);
                 
                    thetaon = problem->DEIMObject->onlineCoeffs(T);
                    Eigen::VectorXd c = problem->DEIMObject->MatrixOnline*thetaon;
                    volVectorField T2("Tonline", Foam2Eigen::Eigen2field(T, c));
                    Info << "###########################################################" << endl;
        #include "createUfIfPresent.H"
                    problem->turbulence->validate();
                    /// Exporting initial solution
                    ITHACAstream::exportSolution(T2, name(counter), folder);
                    ITHACAstream::exportSolution(U, name(counter), folder);
                    ITHACAstream::exportSolution(p, name(counter ), folder);
                    
                    ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
                    std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
                    Ured.append(U.clone());
                    Tred.append(T2.clone());
                    Pred.append(p.clone());
                    scalar  cumulativeContErr = 0.0;
                    counter++;
                    nextWrite += writeEvery;
          
                    while(runTime.run()){

                        runTime.setEndTime(finalTime);
                        runTime++;
                        //p.storePrevIter();
                        Info<< "Time = " << runTime.timeName() << nl << endl;
                        // --- Pressure-velocity PIMPLE corrector loop
                        while (pimple.loop())
                        {
                                 //#include "UEqn.H"

                            fvVectorMatrix UEqn1
                            (
                               fvm::ddt(U) + problem->turbulence->divDevReff(U) 
                               == 
                               T2 + fvOptions(U)  // T2 = deim(div(phi, U)); div(phi, U) evaluated at deim points
                            );
                            
                            UEqn1.relax();
                            fvOptions.constrain(UEqn1);
                            List<Eigen::MatrixXd> RedLinSysU;
                            
                            if (pimple.momentumPredictor())
                            {
                                //solve(UEqn1 == -fvc::grad(p));
                                //UEqn1.solve();
                                RedLinSysU = problem->Umodes.project(UEqn1, NmodesU);
                                volVectorField gradpfull = -fvc::grad(p);
                                Eigen::MatrixXd projGrad =  problem->Umodes.project(gradpfull, NmodesU);
                                RedLinSysU[1] = RedLinSysU[1] + projGrad;
                                a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual);
                                problem->Umodes.reconstruct(U, a, "U");                                
                                fvOptions.correct(U);
                            }

                            // --- Pressure corrector loop
                            while (pimple.correct())
                            {
                                //#include "pEqn.H"
                                volScalarField rAU(1.0 / UEqn1.A());
                                volVectorField HbyA(constrainHbyA(rAU * UEqn1.H(), U, p)); //p
                                surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA)   + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi));

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
                                //bold = b;

                                // Non-orthogonal pressure corrector loop
                                while (pimple.correctNonOrthogonal())
                                {
                                    fvScalarMatrix pEqn
                                    (
                                        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
                                    );

                                    pEqn.setReference(pRefCell, pRefValue);
                                    //pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
                                    RedLinSysP = problem->Pmodes.project(pEqn, NmodesP);
                                    b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
                                    problem->Pmodes.reconstruct(p, b, "p");

                                    if (pimple.finalNonOrthogonalIter())
                                    {
                                        phi = phiHbyA - pEqn.flux();
                                    }
                                }
                                // Explicitly relax pressure for momentum corrector
                                //p.relax();
                                U = HbyA - rAtU * fvc::grad(p); //p
                                U.correctBoundaryConditions();
                                fvOptions.correct(U);
                            } // End of the pimple correct
                            
                            if (pimple.turbCorr())
                            {
                                laminarTransport.correct();
                                problem->turbulence->correct();
                            }  
                        }// end pimple loop
                        problem->Umodes.reconstruct(U, a, "U");
                        problem->Pmodes.reconstruct(p, b, "p");
                        //T = -fvc::div(phi, U); // -fvc::grad(p);
                        T = -fvc::div(phi, U); // - fvc::grad(p);  
                        thetaon = problem->DEIMObject->onlineCoeffs(T);
                        c = problem->DEIMObject->MatrixOnline*thetaon;
                        T2 =  Foam2Eigen::Eigen2field(T, c);
                        if(checkWrite(runTime))
                        {
                            ITHACAstream::exportSolution(U, name(counter), folder);
                            ITHACAstream::exportSolution(T2, name(counter), folder);
                            ITHACAstream::exportSolution(p, name(counter), folder);
                            ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
                            std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
                            Ured.append(U.clone());
                            Pred.append(p.clone());
                            Tred.append(T2.clone());
                            counter++;
                            nextWrite += writeEvery;
                        }        
                       
                       
                    }// end runTime loop

                }

};    

        
int main(int argc, char* argv[])
{
    
    // Construct the tutorial object
    //autoPtr<FSIDEIM> fsideim = autoPtr<FSIDEIM>(new FSIDEIM(argc, argv));
    FSIDEIM fsideim(argc, argv);
    ITHACAparameters* para = ITHACAparameters::getInstance(fsideim.meshPtr(), fsideim._runTime());
    fsideim.startTime = 0.0;
    fsideim.finalTime = 5;
    fsideim.timeStep = 0.0025;
    fsideim.writeEvery = 0.025;
    //Perform the offline solve
    fsideim.offlineSolve();
    ITHACAPOD::getModes(fsideim.Ufield,   fsideim.Umodes,   fsideim._U().name(),0, 0, 0, fsideim.NmodesU);
    ITHACAPOD::getModes(fsideim.Pfield,   fsideim.Pmodes,   fsideim._p().name(),   0, 0, 0, fsideim.NmodesP);
    //ITHACAPOD::getModes(fsideim.Tfield, fsideim.Tmodes, fsideim._T().name(), 0, 0, 0, fsideim.NmodesS);

    // Compute the offline part of the DEIM procedure
    fsideim.PODDEIM(fsideim.NmodesS);
    // restart the simulation
    fsideim.restart();

    //return 0;
    ReducedFSIDEIM test(fsideim);
    test.startTime = 0.0;
    test.finalTime = 5;
    test.timeStep = 0.0025; //0.01;
    test.writeEvery = 0.025;
    test.onlineSolve(fsideim.NmodesU, fsideim.NmodesP,  fsideim.NmodesS);
    std::cout << "################ Error of U #####################" << std::endl;
    Eigen::MatrixXd errL2u = ITHACAutilities::errorL2Rel(fsideim.Ufield, test.Ured);
    std::cout << "################  Error of div(phi, U) #####################" << std::endl;
    Eigen::MatrixXd errL2c = ITHACAutilities::errorL2Rel(fsideim.Tfield, test.Tred);
std::cout << "############### Error p  ######################" << std::endl;
    Eigen::MatrixXd errL2p = ITHACAutilities::errorL2Rel(fsideim.Pfield, test.Pred);
    ITHACAstream::exportMatrix(errL2u, "errL2u", "python","./ITHACAoutput/ErrorsL2/");
    ITHACAstream::exportMatrix(errL2c, "errL2c", "python","./ITHACAoutput/ErrorsL2/");
    ITHACAstream::exportMatrix(errL2p, "errL2p", "python","./ITHACAoutput/ErrorsL2/");
return 0;

//     return 0;
}
