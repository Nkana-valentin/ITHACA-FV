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
//#include <tuple>
//#include <iostream>
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
#define _USE_MATH_DEFINES
#include <cmath>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include <tuple>



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
            theta = P.transpose()*ftoeigen;
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

//*******************************************************************************
//                 FULL ORDER PROBLEM DEFINITION
//*******************************************************************************
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
            NU  = readInt(ITHACAdict->lookup("NU"));
            NS  = readInt(ITHACAdict->lookup("NS")); 
            NP  = readInt(ITHACAdict->lookup("NP")); 
            NR  = readInt(ITHACAdict->lookup("NR")); 

            // second order time derivative
            timeDerivativeSchemeOrder = para->ITHACAdict->lookupOrDefault<word>("timeDerivativeSchemeOrder", "second");

            //std::cout << "####### test of the fsideim ctor ######" << std::endl;
        }

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
        PtrList<volVectorField> Sfield; 
        PtrList<volVectorField> Tfield; 
        PtrList<volVectorField> Rfield; 
        PtrList<volScalarField> Pfield;
        
        autoPtr<IOMRFZoneList> _MRF;
        autoPtr<volScalarField> _p;
        autoPtr<volVectorField> _U;
        autoPtr<volVectorField> _T;
        autoPtr<volVectorField> _S; //to forcing term (fvc::grad(p))
        autoPtr<volVectorField> _R; //list of the residual
        autoPtr<pointVectorField> _pd;
        autoPtr<surfaceVectorField> _Uf; 
        autoPtr<surfaceScalarField> _phi;
        autoPtr<surfaceScalarField> _phi0;
        autoPtr<Foam::dynamicFvMesh> meshPtr;
        bool  correctPhi;
        bool  checkMeshCourantNo;
        bool  moveMeshOuterCorrectors;
        ITHACAparameters* para;
        vectorField point0;
        autoPtr<singlePhaseTransportModel> _laminarTransport;
        autoPtr<incompressible::turbulenceModel> turbulence;
        autoPtr<sixDoFRigidBodyMotionSolver> sDRBMS;
        IOdictionary* dyndict;
        vVf* DEIMObject;

        volVectorModes Smodes;
        volScalarModes Pmodes;
        // List of Modes for U
        volVectorModes Umodes;
        // List of modes for the residual
         volVectorModes Rmodes;
        // List of Modes for Phi
        surfaceScalarModes Phimodes;
        // Number of Modes to approximate U, S and p
        int NU, NS, NP, NR;
         // Reduced source terms
        Eigen::MatrixXd ReducedVectors;
        // POD modes in Eigen Format
        Eigen::MatrixXd ModesUEig;
        label counter = 0;
        // second order time derivative
        word timeDerivativeSchemeOrder;

        /// Offline problem 
        void offlineSolve(fileName folder = "./ITHACAoutput/Offline") 
        {

            Time& runTime = _runTime();
            surfaceScalarField& phi = _phi();
            dynamicFvMesh& mesh = meshPtr();
            fv::options& fvOptions = _fvOptions();
            pimpleControl& pimple = _pimple();
            volScalarField& p = _p();
            volVectorField& U = _U();
            volVectorField& S = _S();
            volVectorField& T = _T();
            volVectorField& R = _R();
            surfaceVectorField&  Uf = _Uf(); 
            IOMRFZoneList& MRF = _MRF();
            singlePhaseTransportModel& laminarTransport = _laminarTransport();
            turbulence = autoPtr<incompressible::turbulenceModel>(incompressible::turbulenceModel::New(U, phi, laminarTransport));
            label pRefCell = 0;
            scalar pRefValue = 0.0;
            volScalarField nueff = laminarTransport.nu();
            instantList Times = runTime.times();
            runTime.setEndTime(finalTime);
            runTime.setTime(Times[1], 1);
            runTime.setDeltaT(timeStep);
            nextWrite = startTime; 
            label counter = 1;
//#include "createUfIfPresent.H"
            turbulence->validate();
            /// Exporting initial solution
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(R, name(counter), folder);
            ITHACAstream::exportSolution(T, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter ), folder);
            
            ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
            std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
            Ufield.append(U.clone());
            //Sfield.append(S.clone());
            Tfield.append(T.clone());
            Rfield.append(R.clone());
            Pfield.append(p.clone());
            //Phifield.append(phi.clone());
            scalar  cumulativeContErr = 0.0;
            counter++;
            nextWrite += writeEvery;
  
            while(runTime.run()){

                runTime.setEndTime(finalTime);
                runTime++;
                Info<< "Time = " << runTime.timeName() << nl << endl;
       
                while (pimple.loop())
                {
                    // if (pimple.firstIter() || moveMeshOuterCorrectors)
                    // {
                    //     // Do any mesh changes
                    //     mesh.controlledUpdate();

                    //     if (mesh.changing())
                    //     {
                    //         MRF.update();

                    //         if (correctPhi)
                    //         {
                    //             // Calculate absolute flux
                    //             // from the mapped surface velocity
                    //             phi = mesh.Sf() & Uf;

                    //             #include "correctPhi.H"

                    //             // Make the flux relative to the mesh motion
                    //             fvc::makeRelative(phi, U);
                    //         }

                    //         if (checkMeshCourantNo)
                    //         {
                    //             #include "meshCourantNo.H"
                    //         }
                    //     }
                    // }

                    #include "UEqn.H"

                    // --- Pressure corrector loop
                    while (pimple.correct())
                    {
                        #include "pEqn1.H"
                    }

                    if (pimple.turbCorr())
                    {
                        laminarTransport.correct();
                        turbulence->correct();
                    }
                }

                // Eval the non-linear term  
                T = -fvc::div(phi,U);       
                if(checkWrite(runTime))
                {
                    ITHACAstream::exportSolution(U, name(counter), folder);
                    ITHACAstream::exportSolution(T, name(counter), folder);
                    ITHACAstream::exportSolution(R, name(counter), folder);
                    //ITHACAstream::exportSolution(S, name(counter), folder); //snapshots of the forcing term
                    ITHACAstream::exportSolution(p, name(counter), folder);
                    ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
                    std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
                    Ufield.append(U.clone());
                    Pfield.append(p.clone());
                    Tfield.append(T.clone());
                    Rfield.append(R.clone());
                    counter++;
                    nextWrite += writeEvery;
                }        
               
               
            }// end runTime loop

        }
  
        void PODDEIM(int NmodesS)
        {
           
            /// Create DEIM object with given number of basis functions
             dynamicFvMesh& mesh = meshPtr();
             DEIMObject = new vVf(Rfield, NR, "residualterm", _R().name());
             DEIMObject->subfield = autoPtr<volVectorField>(new volVectorField(DEIMObject->generateSubmesh(1, mesh, _R())));
           
        }
       

        bool checkWrite(Time& timeObject)
        {
            scalar diffnow = mag(nextWrite - atof(timeObject.timeName().c_str()));
            scalar diffnext = mag(nextWrite - atof(timeObject.timeName().c_str()) -  timeObject.deltaTValue());

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
                 _T.clear();
                 _p.clear();
                _U.clear();
                _phi.clear();
                //_Uf.clear();
                _pd.clear();
               
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
/// Structure for the resolution of the nonlinear problem
struct newtonObject: public newton_argument<double>
{
    public:
        newtonObject() {}
        newtonObject(int Nx, int Ny,
                     FSIDEIM& problem): newton_argument<double>(Nx, Ny), problem(&problem)
        {}
        // function which computes the residual
        int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
        // function that computes the Jacobian
        int df(const Eigen::VectorXd& x,  Eigen::MatrixXd& fjac) const;

        // pointer to the FOM problem
        FSIDEIM* problem;
        // Time step in the reduced problem
        scalar dt;
        // Reduced solution at the previous time step
        Eigen::VectorXd y_old;

        // second order time derivative
        Eigen::VectorXd y_oldold;

        // Total residual Matrix
        Eigen::MatrixXd R;
        Eigen::MatrixXd T; // from the temporal term
        Eigen::MatrixXd C; //convective term
        Eigen::MatrixXd D; // diffusive term
        //Reduced matrix from the pressure equation
        Eigen::MatrixXd Pr;

        // Reduced source term from pressure equation
        Eigen::VectorXd r;
    
        // total reduced source term
        Eigen::VectorXd s;

};
// function that computes the reduced residual
int newtonObject::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
{
    // x is the tentative value, fvec is the residual
    // Time derivative of reduced coefficients
    Eigen::VectorXd a_dot(problem->NU);
    a_dot = (x - y_old) / dt;

    // second order time derivative
    if (problem->timeDerivativeSchemeOrder == "first")
     {
        a_dot = (x - y_old) / dt;
     }
     else
     {
         a_dot = (1.5 * x - 2 * y_old + 0.5 * y_oldold) / dt;
     }

     fvec = T*a_dot + C*x + D*x - s; //total residual
     //fvec =     R*a_dot - s; //total residual
     return 0;
}

// function that computes the Jacobian of the NL problem, it is done by numerical differentiation
int newtonObject::df(const Eigen::VectorXd& x,Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newtonObject> numDiff(*this);
    numDiff.df(x, fjac);
    std::cout << "============================ inside df =================================" << std::endl;
    return 0;
}



/// class for the reduced problem
class ReducedFSIDEIM : public reductionProblem
{
        public:
                explicit ReducedFSIDEIM(FSIDEIM& FoamProblem) : problem(&FoamProblem)
                {

           std::cout << "============================ inside reduced deim =================================" << std::endl;
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
                     a = Eigen::VectorXd::Zero(problem->NU);
                     b = Eigen::VectorXd::Zero(problem->NP);
                     newton_object = newtonObject(problem->NU, problem->NU, FoamProblem);
                     newton_objectP = newtonObject(problem->NP, problem->NP, FoamProblem);

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
                // Reduced vectors
                Eigen::VectorXd a;
                Eigen::VectorXd b;
                // Reduced residual vectors
                Eigen::VectorXd ures;
                Eigen::VectorXd pres;  

                /// Function object to call the non linear solver sup approach
                newtonObject newton_object;
                newtonObject newton_objectP;

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
                    volVectorField& R = problem->_R();
                    IOMRFZoneList& MRF = problem->_MRF();
                    surfaceVectorField&  Uf = problem->_Uf(); 
                    bool  correctPhi = problem->correctPhi;
                    bool  checkMeshCourantNo = problem->checkMeshCourantNo;
                    bool  moveMeshOuterCorrectors = problem->moveMeshOuterCorrectors;
                    sixDoFRigidBodyMotionSolver& sDRBMS = problem->sDRBMS();
                    singlePhaseTransportModel& laminarTransport = problem->_laminarTransport();
                    volScalarField nueff = laminarTransport.nu();
                    label pRefCell = 0;
                    scalar pRefValue = 0.0;
                  
                    instantList Times = runTime.times();
                    runTime.setEndTime(finalTime);
                    runTime.setTime(Times[1], 1);
                    runTime.setDeltaT(timeStep);
                    nextWrite = startTime; 
                    label counter = 1;
                    Eigen::VectorXd a0 = Eigen::VectorXd::Zero(NmodesU);
                    Eigen::VectorXd b = Eigen::VectorXd::Zero(NmodesU);
                   
                    a0 = ITHACAutilities::getCoeffs(U, problem->Umodes, NmodesU, true);
                    std::cout << a0 << std::endl;
                    std::cout << "============================ a0 =============================="<< std::endl;
                      //Evaluate the new forcing term using DEIM
                    //Eigen::VectorXd  thetaon = problem->DEIMObject->onlineCoeffs(S);
                    Eigen::VectorXd  thetaon = problem->DEIMObject->onlineCoeffs(R);
                        // newton_object.s = problem->DEIMObject->MatrixOnline*thetaon;
                         // newton_object.s = problem->ReducedVectors*thetaon;
                        // Interpolate the forcing term
                        Eigen::VectorXd  r = problem->DEIMObject->MatrixOnline*thetaon; 
                        volVectorField R2("Sonline", Foam2Eigen::Eigen2field(R, r)); //declare the field outside the loop
                    //b = ITHACAutilities::getCoeffs(p, problem->Pmodes, NmodesP, true);
//problem->restart();
                    // set initial condition to zero with a0
                    a(a0);
                   // Eigen::VectorXd a = Eigen::VectorXd::Zero(NmodesU);
                    newton_object.y_old = a;
                    // second order time derivative
                    newton_object.y_oldold = a;
                    newton_object.dt = timeStep;
                    Color::Modifier red(Color::FG_RED);
                    Color::Modifier green(Color::FG_GREEN);
                    Color::Modifier def(Color::FG_DEFAULT);         
        //#include "createUfIfPresent.H"
                    problem->turbulence->validate();
                    /// Exporting initial solution
                    ITHACAstream::exportSolution(R2, name(counter), folder);
                    ITHACAstream::exportSolution(U, name(counter), folder);
                    ITHACAstream::exportSolution(p, name(counter ), folder);
                    
                    ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
                    std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
                    Ured.append(U.clone());
                    //Tred.append(T2.clone());
                    Pred.append(p.clone());
                    scalar  cumulativeContErr = 0.0;
                    counter++;
                    nextWrite += writeEvery;
                    int i=0;
                    label time_rom=0;
          
                    while(runTime.run()){

                        runTime.setEndTime(finalTime);
                        runTime++;
                        //p.storePrevIter();
                        Info<< "Time = " << runTime.timeName() << nl << endl;
                        // --- Pressure-velocity PIMPLE corrector loop
                        while (pimple.loop())
                        {
                            // if (pimple.firstIter() || moveMeshOuterCorrectors)
                            // {
                            //     // Do any mesh changes
                            //     mesh.controlledUpdate();

                            //     if (mesh.changing())
                            //     {
                            //         MRF.update();

                            //         if (correctPhi)
                            //         {
                            //             // Calculate absolute flux
                            //             // from the mapped surface velocity
                            //             phi = mesh.Sf() & Uf;

                            //             #include "correctPhi.H"

                            //             // Make the flux relative to the mesh motion
                            //             fvc::makeRelative(phi, U);
                            //         }

                            //         if (checkMeshCourantNo)
                            //         {
                            //             #include "meshCourantNo.H"
                            //         }
                            //     }
                            // }
                            //#include "UEqn.H"
                            fvVectorMatrix UEqn
                            (
                               fvm::ddt(U) 
                               + fvm::div(phi, U)
                               +  problem->turbulence->divDevReff(U)                            
                               //== fvOptions(U) 
                            );
                            fvVectorMatrix UEqn1( fvm::ddt(U) );
                            fvVectorMatrix UEqn2( fvm::div(phi, U));
                            fvVectorMatrix UEqn3( problem->turbulence->divDevReff(U) );
                            volVectorField sourcefull = -fvc::grad(p); // + fvOptions(U);?????
                            //RedLinSysU[1] = RedLinSysU[1] + projGrad;

                            vectorField resF(UEqn.residual());

                            for (label i = 0; i < R.internalField().size(); i++)
                            {
                                R.ref()[i] = resF[i];
                            }
                                                        
                            UEqn.relax();
                            fvOptions.constrain(UEqn);
                            //UEqn1.relax();
                            //fvOptions.constrain(UEqn1);
                            UEqn2.relax();
                            fvOptions.constrain(UEqn2);
                            UEqn3.relax();
                            fvOptions.constrain(UEqn3);
                            // List of matrices for reduced systems.
                            List<Eigen::MatrixXd> projfvMat;///???????
                            List<Eigen::MatrixXd> projfvMat1;
                            List<Eigen::MatrixXd> projfvMat2;
                            List<Eigen::MatrixXd> projfvMat3;
                            Eigen::MatrixXd rhs; //rhs of UEqn 

                            if (pimple.momentumPredictor())
                            {
                           
                                projfvMat  = problem->Umodes.project(UEqn,   NmodesU);
                                projfvMat1 = problem->Umodes.project(UEqn1, NmodesU);
                                projfvMat2 = problem->Umodes.project(UEqn2, NmodesU);
                                projfvMat3 = problem->Umodes.project(UEqn3, NmodesU);
                                rhs = problem->Umodes.project(sourcefull,   NmodesU);

                                newton_object.R = projfvMat[0];
                                newton_object.T = projfvMat1[0];
                                newton_object.C = projfvMat2[0];
                                newton_object.D = projfvMat3[0];
                                // Total reduced rhs
                                newton_object.s = projfvMat1[1] + projfvMat2[1] + projfvMat3[1] + rhs;
                                //newton_object.s = projfvMat[1] + rhs;
                                std::cout << "============================ line 661 =============================="<< std::endl;
                                Eigen::HybridNonLinearSolver<newtonObject> hnls(newton_object);
                                // create empty vector as a copy of a
                                Eigen::VectorXd res(a);
                                std::cout << "============================ line 668 =============================="<< std::endl;
                                auto t1 = std::chrono::high_resolution_clock::now();
                                // solve the reduced problem
                                hnls.solve(a);
                                std::cout << "============================ line 673 =============================="<< std::endl;

                                auto t2 = std::chrono::high_resolution_clock::now();
                                auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
                                time_rom +=time_span.count();
                                // compute the residual to check accuracy
                                newton_object.operator()(a, res);
                            std::cout << "============================ line 680 =============================="<< std::endl;

                                // second order time derivative
                                newton_object.y_oldold = newton_object.y_old;
                                // set y_old=a
                                newton_object.y_old = a;
                                // Output some statics on the nonlinear solver
                                if (res.norm() < 1e-5)
                                {
                                    std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                                              hnls.iter << " iterations " << def << std::endl << std::endl;
                                }

                                else
                                {
                                    std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                                              hnls.iter << " iterations " << def << std::endl << std::endl;
                                }
                                // Reconstruct solution  
                                volVectorField U("U", problem->Umodes[0]); // Initialization of U before the reconstruction
                                problem->Umodes.reconstruct(U, a, "U");                             
                                fvOptions.correct(U);
                                //return ;
                            }

                            // --- Pressure corrector loop
                            while (pimple.correct())
                            {
                                //#include "pEqn.H"
                                volScalarField rAU(1.0 / UEqn.A());
                                volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p)); //p
                                surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA)+fvc::interpolate(rAU)*fvc::ddtCorr(U, Uf));

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
                                    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
                                    //RedLinSysP = problem->Pmodes.project(pEqn, NmodesP);
                                    //b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
                                    //problem->Pmodes.reconstruct(p, b, "p");

                                    if (pimple.finalNonOrthogonalIter())
                                    {
                                        phi = phiHbyA - pEqn.flux();
                                    }
                                }
                                // Explicitly relax pressure for momentum corrector
                                //p.relax();

                                U = HbyA - rAtU*fvc::grad(p);
                                U.correctBoundaryConditions();
                                fvOptions.correct(U);
                                {
                                    Uf = fvc::interpolate(U);
                                    surfaceVectorField n(mesh.Sf()/mesh.magSf());
                                    Uf += n*(phi/mesh.magSf() - (n & Uf));
                                }

                                // Make the fluxes relative to the mesh motion
                                fvc::makeRelative(phi, U);
                            } // End of the pimple correct
                            
                            if (pimple.turbCorr())
                            {
                                laminarTransport.correct();
                                problem->turbulence->correct();
                            }  
                        }// end pimple loop

                        if(checkWrite(runTime))
                        {
                            ITHACAstream::exportSolution(U, name(counter), folder);
                            //ITHACAstream::exportSolution(T2, name(counter), folder);
                            ITHACAstream::exportSolution(p, name(counter), folder);
                            ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
                            std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
                            Ured.append(U.clone());
                            Pred.append(p.clone());
                            //Tred.append(T2.clone());
                            counter++;
                            nextWrite += writeEvery;
                        }        
                       
                       
                    }// end runTime loop

                }
                bool checkWrite(Time& timeObject)
                {
                    scalar diffnow = mag(nextWrite - atof(timeObject.timeName().c_str()));
                    scalar diffnext = mag(nextWrite - atof(timeObject.timeName().c_str()) -  timeObject.deltaTValue());

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
       
int main(int argc, char* argv[])
{
    
    // Construct the tutorial object
    FSIDEIM fsideim(argc, argv);
    ITHACAparameters* para = ITHACAparameters::getInstance(fsideim.meshPtr(), fsideim._runTime());
    fsideim.startTime = 0.0;
    fsideim.finalTime = 1;
    fsideim.timeStep = 0.001;
    fsideim.writeEvery = 0.01;
    //Perform the offline solve
    fsideim.offlineSolve();
    ITHACAPOD::getModes(fsideim.Ufield,   fsideim.Umodes,   fsideim._U().name(),   0, 0, 0,    fsideim.NU);
    ITHACAPOD::getModes(fsideim.Pfield,   fsideim.Pmodes,   fsideim._p().name(),   0, 0, 0,    fsideim.NP);
    //ITHACAPOD::getModes(fsideim.Rfield,   fsideim.Rmodes,   fsideim._R().name(),   0, 0, 0,    fsideim.NR);
    //ITHACAPOD::getModes(fsideim.Sfield, fsideim.Smodes, fsideim._S().name(), 0, 0, 0, fsideim.NS);
    //return 0;
    // Compute the offline part of the DEIM procedure
    fsideim.PODDEIM(fsideim.NR);
std::cout << "################ "<< "End of the Offline part "<< " ##################" << std::endl;
//return 0;

    fsideim.restart();


    ReducedFSIDEIM test(fsideim);
    test.startTime =  fsideim.startTime;
    test.finalTime =  fsideim.finalTime;
    test.timeStep  =  fsideim.timeStep ; //0.01;
    test.writeEvery = fsideim.writeEvery;
    test.onlineSolve(fsideim.NU, fsideim.NP,  fsideim.NS);

   // Eigen::MatrixXd errL2s = ITHACAutilities::errorL2Rel(fsideim.Tfield, test.Tred);
   // std::cout << "##################################" << std::endl;
   Eigen::MatrixXd errL2u = ITHACAutilities::errorL2Rel(fsideim.Ufield, test.Ured);

   //  std::cout << "##################################" << std::endl;
   //  Eigen::MatrixXd errL2p = ITHACAutilities::errorL2Rel(fsideim.Pfield, test.Pred);
return 0;

//     return 0;
}


















//void PODDEIM(int NmodesU, int NmodesS)
       // {
           //      NU = NmodesU;
           //      NS = NmodesS;
           //      Mr.resize(NU, NU);
           //      Br.resize(NU, NU);
           //      Cr.resize(NU, NS);
           //      Sr.resize(NU, NS);
           //      Hr.resize(NU, 1);
           //      volScalarField nueff = _laminarTransport().nu();

           //      for (int i = 0; i < NU; i++)
           //      {
           //          for (int j = 0; j < NU; j++)
           //          {
           //              Mr(i, j) =  fvc::domainIntegrate(Umodes[i] & Umodes[j]).value(); 
           //              Br(i, j) = -fvc::domainIntegrate(Umodes[i] & fvc::laplacian(nueff, Umodes[j])).value() 
           //                         -fvc::domainIntegrate(Umodes[i] & fvc::div(nueff * dev2(T(fvc::grad(Umodes[j])))) ).value();
           //              //Might be computed online?????
           //              //Cr(i, j) = fvc::domainIntegrate(Umodes[i] & fvc::div(_phi(), Umodes[j])).value();
           //          }
           //      }
                
           //   std::cout << Mr << std::endl;
           // std::cout << "============================ Mr =================================" << std::endl;
           //      std::cout << Br << std::endl;

           // std::cout << "============================ Br =================================" << std::endl;

           // for (label i = 0; i < NU; i++)
           // { 
           //      Hr(i) = -fvc::domainIntegrate(Umodes[i] & fvc::grad(_p()) ).value();   
           //  }

           //   std::cout << Hr << std::endl;

           // std::cout << "============================  Hr =================================" << std::endl;
                
            /// Create DEIM object with given number of basis functions
             // dynamicFvMesh& mesh = meshPtr();
             // DEIMObject = new vVf(Tfield, NS, "convectiveterm", _T().name());
             // DEIMObject->subfield = autoPtr<volVectorField>(new volVectorField(DEIMObject->generateSubmesh(1, mesh, _T())));
             // for (label i = 0; i < NU; i++)
             // { 
             //       for (label j = 0; j < NS; j++)
             //       {
             //            //Sr(i, j) = -fvc::domainIntegrate(Umodes[i] & DEIMObject->modes[j] ).value();
             //            Cr(i, j) = fvc::domainIntegrate(Umodes[i] & DEIMObject->modes[j] ).value();///-???
             //       }    
             // }

             //std::cout << Cr << std::endl;

           //std::cout << "============================  Cr =================================" << std::endl;
               


             //ModesUEig = Foam2Eigen::PtrList2Eigen(Umodes);
             //ModesUEig.conservativeResize(ModesUEig.rows(), NU);
 //            auto vol = ITHACAutilities::getMassMatrixFV(_S());

 //            Info << vol << endl;
 //             std::cout << "============================ vol =================================" << std::endl;
 //             ReducedVectors = ((DEIMObject->P).transpose()*(DEIMObject->U)).fullPivLu().inverse();
 // std::cout << ReducedVectors << std::endl;
 //              std::cout << "============================  ReducedVectors =================================" << std::endl;
           
        //}
