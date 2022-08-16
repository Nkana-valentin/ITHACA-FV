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


class DEIM_functionU : public DEIM<fvVectorMatrix>
{
    public:
        using DEIM::DEIM;
        autoPtr<volVectorField> fieldsA;
        autoPtr<volVectorField> fieldsB;
        autoPtr<surfaceScalarField> phi;

        Eigen::MatrixXd onlineCoeffsAU(fvVectorMatrix& Ueqn)
        {
            std::cout << "#################" << magicPointsAcol().size() << "####################" << std::endl; 
            Eigen::MatrixXd theta(magicPointsAcol().size(), 1);
            Eigen::SparseMatrix<double> Mr;
            Eigen::VectorXd br;
            Foam2Eigen::fvMatrix2Eigen(Ueqn, Mr, br);
            for (int i = 0; i < magicPointsAcol().size(); i++)
            {
                int ind_row = magicPointsArow()[i] + xyz_Arow()[i] * fieldsA().size();
                int ind_col = magicPointsAcol()[i] + xyz_Acol()[i] * fieldsA().size();
                theta(i) = Mr.coeffRef(ind_row, ind_col);
            }
            std::cout << "#################" << "End onlineCoeffsAU method" << "####################" << std::endl; 
            return theta;
            
        }

        Eigen::MatrixXd onlineCoeffsBU(fvVectorMatrix& Ueqn)
        {
            Eigen::MatrixXd theta(magicPointsB().size(), 1);
            //Ueqn.relax(); //can be cancelled assuming that we relax the ueqn argument
            Eigen::SparseMatrix<double> Mr;
            Eigen::VectorXd br;
            Foam2Eigen::fvMatrix2Eigen(Ueqn, Mr, br);
            for (int i = 0; i < magicPointsB().size(); i++)
            {
                int ind_row = magicPointsB()[i] + xyz_B()[i] *  fieldsB().size();
                theta(i) = br(ind_row);
            }
            std::cout << "#################" << "End onlineCoeffsBU method" << "####################" << std::endl; 
            return theta;
            
        }

        
};
class DEIM_functionP : public DEIM<fvScalarMatrix>
{
    public:
        using DEIM::DEIM;
        autoPtr<volScalarField> fieldsA;
        autoPtr<volScalarField> fieldsB;
        autoPtr<surfaceScalarField> phi;
        Eigen::MatrixXd onlineCoeffsAP(fvScalarMatrix& peqn)
        {
            std::cout << "#################" << "start onlineCoeffsAP method" << "####################" << std::endl; 
            Eigen::MatrixXd theta(magicPointsAcol().size(), 1);
            Eigen::SparseMatrix<double> Mr;
            Eigen::VectorXd br;
            Foam2Eigen::fvMatrix2Eigen(peqn, Mr, br);
            for (int i = 0; i < magicPointsAcol().size(); i++)
            {
                int ind_row = magicPointsArow()[i] + xyz_Arow()[i] * fieldsA().size();
                int ind_col = magicPointsAcol()[i] + xyz_Acol()[i] * fieldsA().size();
                theta(i) = Mr.coeffRef(ind_row, ind_col);
            }
            std::cout << "#################" << "End onlineCoeffsAP method" << "####################" << std::endl; 
            return theta;

         
        }
        Eigen::MatrixXd onlineCoeffsBP(fvScalarMatrix& peqn)
        {
                Eigen::MatrixXd theta(magicPointsB().size(), 1);
                Eigen::SparseMatrix<double> Mr;
                Eigen::VectorXd br;
                Foam2Eigen::fvMatrix2Eigen(peqn, Mr, br);
                for (int i = 0; i <  magicPointsB().size(); i++)
                {
                    int ind_row = magicPointsB()[i] + xyz_B()[i] * fieldsB().size();
                    theta(i) = br(ind_row);
                }
                std::cout << "#################" << "End onlineCoeffsBP method" << "####################" << std::endl; 

                return theta;
        }
        
};
class DEIM_Ufield : public DEIM<volVectorField>
{
    public:
        using DEIM::DEIM;
        autoPtr<volVectorField> fieldsA;
        autoPtr<volVectorField> fieldsB;
        autoPtr<surfaceScalarField> phi;


        
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
            NmodesP  = readInt(ITHACAdict->lookup("NmodesP"));
            NmodesDEIMAU = readInt(ITHACAdict->lookup("NmodesDEIMAU")); 
            NmodesDEIMbU = readInt(ITHACAdict->lookup("NmodesDEIMbU")); 
            NmodesDEIMAP = readInt(ITHACAdict->lookup("NmodesDEIMAP"));
            NmodesDEIMbP = readInt(ITHACAdict->lookup("NmodesDEIMbP"));


            std::cout << "####### test of the fsideim ctor ######" << std::endl;
        }
        //FSIDEIM() {};
        ~FSIDEIM() {};

        /// Fields To Perform
        /// Velocity field
        // volVectorField& U;
        //  /// Pressure field
        // volScalarField& p;
        // /// flux field
        // surfaceScalarField& phi;
        // /// pointDisplacement field
        //pointVectorField& pd;
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
        PtrList<volScalarField> Pfield;
        PtrList<surfaceScalarField> Phifield;
        
        // Pointer to the residual solution
        PtrList<fvVectorMatrix> resField;
        autoPtr<IOMRFZoneList> _MRF;
        autoPtr<volScalarField> _p;
        autoPtr<volVectorField> _U;
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

        PtrList<fvScalarMatrix> pEqnlist;
        PtrList<fvVectorMatrix> UEqnlist;
        PtrList<fvVectorMatrix> DivList;
        PtrList<surfaceScalarField> Philist;

        Eigen::MatrixXd ModesUEig;
        Eigen::MatrixXd ModesPEig;
        Eigen::MatrixXd ModesPhiEig;
     
        // autoPtr<DEIM_functionU> _DEIMmatriceU; 
        // autoPtr<DEIM_functionP> _DEIMmatriceP; 
        DEIM_functionU* DEIMmatriceU; 
        DEIM_functionP* DEIMmatriceP; 
        autoPtr<DEIM_Ufield> DEIMfieldU;

        std::vector<Eigen::MatrixXd> ReducedMatricesA_U;
        std::vector<Eigen::MatrixXd> ReducedVectorsB_U;
        std::vector<Eigen::MatrixXd> ReducedMatricesA_p;
        std::vector<Eigen::MatrixXd> ReducedVectorsB_p;

        PtrList<volScalarField> pFieldsA_U;
        PtrList<volVectorField> UFieldsA_U;
        PtrList<surfaceScalarField> phiFieldsA_U;
        PtrList<volScalarField> pFieldsB_U;
        PtrList<volVectorField> UFieldsB_U;
        PtrList<surfaceScalarField> phiFieldsB_U;
        PtrList<volScalarField> pFieldsA_p;
        PtrList<volVectorField> UFieldsA_p;
        PtrList<surfaceScalarField> phiFieldsA_p;
        PtrList<volScalarField> pFieldsB_p;
        PtrList<volVectorField> UFieldsB_p;
        PtrList<surfaceScalarField> phiFieldsB_p;
        // List of Modes for P
        volScalarModes Pmodes;
        // List of Modes for U
        volVectorModes Umodes;
        // List of Modes for Phi
        surfaceScalarModes Phimodes;

        int NmodesU;
        int NmodesP;
        int NmodesDEIMAU; 
        int NmodesDEIMbU;
        int NmodesDEIMAP;
        int NmodesDEIMbP;  
        Eigen::MatrixXd onlineCoeffsAP(fvScalarMatrix& peqn)
        {
            std::cout << "#################" << "start onlineCoeffsAP method" << "####################" << std::endl; 
            Eigen::MatrixXd theta(DEIMmatriceP->magicPointsAcol().size(), 1);
            Eigen::SparseMatrix<double> Mr;
            Eigen::VectorXd br;
            Foam2Eigen::fvMatrix2Eigen(peqn, Mr, br);
            for (int i = 0; i < DEIMmatriceP->magicPointsAcol().size(); i++)
            {
                int ind_row = DEIMmatriceP->magicPointsArow()[i] + DEIMmatriceP->xyz_Arow()[i] * _p().size();
                int ind_col = DEIMmatriceP->magicPointsAcol()[i] + DEIMmatriceP->xyz_Acol()[i] * _p().size();
                theta(i) = Mr.coeffRef(ind_row, ind_col);
            }
            std::cout << "#################" << "End onlineCoeffsAP method" << "####################" << std::endl; 
            return theta;

         
        }
        Eigen::MatrixXd onlineCoeffsBP(fvScalarMatrix& peqn)
        {
                Eigen::MatrixXd theta(DEIMmatriceP->magicPointsB().size(), 1);
                Eigen::SparseMatrix<double> Mr;
                Eigen::VectorXd br;
                Foam2Eigen::fvMatrix2Eigen(peqn, Mr, br);
                for (int i = 0; i <  DEIMmatriceP->magicPointsB().size(); i++)
                {
                    int ind_row = DEIMmatriceP->magicPointsB()[i] + DEIMmatriceP->xyz_B()[i] * _p().size();
                    theta(i) = br(ind_row);
                }
                std::cout << "#################" << "End onlineCoeffsBP method" << "####################" << std::endl; 

                return theta;
        }

        Eigen::MatrixXd onlineCoeffsAU(fvVectorMatrix& Ueqn)
        {
            Eigen::MatrixXd theta(DEIMmatriceU->magicPointsAcol().size(), 1);
            Eigen::SparseMatrix<double> Mr;
            Eigen::VectorXd br;
            Foam2Eigen::fvMatrix2Eigen(Ueqn, Mr, br);
            for (int i = 0; i < DEIMmatriceU->magicPointsAcol().size(); i++)
            {
                int ind_row = DEIMmatriceU->magicPointsArow()[i] + DEIMmatriceU->xyz_Arow()[i] * _U().size();
                int ind_col = DEIMmatriceU->magicPointsAcol()[i] + DEIMmatriceU->xyz_Acol()[i] * _U().size();
                theta(i) = Mr.coeffRef(ind_row, ind_col);
            }
                std::cout << "#################" << theta.rows() << "####################" << std::endl;
                std::cout << "#################" << theta.cols() << "####################" << std::endl;  
            std::cout << "#################" << "End onlineCoeffsAU method" << "####################" << std::endl; 

            return theta;
            
        }

        Eigen::MatrixXd onlineCoeffsBU(fvVectorMatrix& Ueqn)
        {
            Eigen::MatrixXd theta(DEIMmatriceU->magicPointsB().size(), 1);
            //Ueqn.relax(); //can be cancelled assuming that we relax the ueqn argument
            Eigen::SparseMatrix<double> Mr;
            Eigen::VectorXd br;
            Foam2Eigen::fvMatrix2Eigen(Ueqn, Mr, br);
            for (int i = 0; i < DEIMmatriceU->magicPointsB().size(); i++)
            {
                int ind_row = DEIMmatriceU->magicPointsB()[i] + DEIMmatriceU->xyz_B()[i] * _U().size();
                theta(i) = br(ind_row);
            }
            std::cout << "#################" << "End onlineCoeffsBU method" << "####################" << std::endl; 
            return theta;
            
        }
        
        void offlineSolve(fileName folder = "./ITHACAoutput/Offline") 
        {

            Time& runTime = _runTime();
            surfaceScalarField& phi = _phi();
            dynamicFvMesh& mesh = meshPtr();
            fv::options& fvOptions = _fvOptions();
            pimpleControl& pimple = _pimple();
            volScalarField& p = _p();
            volVectorField& U = _U();
             /// constructing face velocity
            //surfaceVectorField& Uf = *_Uf;
            IOMRFZoneList& MRF = _MRF();
            singlePhaseTransportModel& laminarTransport = _laminarTransport();
            turbulence = autoPtr<incompressible::turbulenceModel>
            (
                     incompressible::turbulenceModel::New(U, phi, laminarTransport)
            );
            volScalarField nuEff = laminarTransport.nu();
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
            ITHACAstream::exportSolution(p, name(counter ), folder);
            //ITHACAstream::exportSolution(phi, name(runTime.timeIndex()), folder);
            ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
            std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
            Ufield.append(U.clone());
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
                }        
                if(checkWrite(runTime))
                {
                    ITHACAstream::exportSolution(U, name(counter), folder);
                    ITHACAstream::exportSolution(p, name(counter), folder);
                    //ITHACAstream::exportSolution(phi, name(runTime.timeIndex()), folder);
                    //ITHACAstream::exportSolution(Flux, name(runTime.timeIndex()), folder);
                    ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
                    std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
                    Ufield.append(U.clone());
                    Pfield.append(p.clone());
                    //Phifield.append(phi.clone());
                    counter++;
                    nextWrite += writeEvery;
                }        
               
               
            }// end runTime loop

        };
  

        void PODDEIM()
        {
           PODDEIM(NmodesU, NmodesP, NmodesDEIMAU, NmodesDEIMbU,NmodesDEIMAP,NmodesDEIMbP);
        }
        void PODDEIM(int NmodesU, int NmodesP,int NmodesDEIMAU,int NmodesDEIMbU,int NmodesDEIMAP,int NmodesDEIMbP)
        {

            dynamicFvMesh& mesh = meshPtr();
            DEIMmatriceU = new DEIM_functionU(UEqnlist, NmodesDEIMAU, NmodesDEIMbU,"U_matrix");
            DEIMmatriceP = new DEIM_functionP(pEqnlist, NmodesDEIMAP, NmodesDEIMbP,"P_matrix");
 
            DEIMmatriceU->fieldsA = autoPtr<volVectorField>(new volVectorField(DEIMmatriceU->generateSubmeshMatrix(2, mesh, _U())));
            DEIMmatriceU->fieldsB = autoPtr<volVectorField>(new volVectorField(DEIMmatriceU->generateSubmeshVector(2, mesh, _U())));
            std::cout << "#####################################" << std::endl; 
            DEIMmatriceP->fieldsA = autoPtr<volScalarField>(new volScalarField(DEIMmatriceP->generateSubmeshMatrix(2, mesh, _p())));
            DEIMmatriceP->fieldsB = autoPtr<volScalarField>(new volScalarField(DEIMmatriceP->generateSubmeshVector(2, mesh, _p())));
            // DEIMmatriceU = autoPtr<DEIM_functionU>(new DEIM_functionU(UEqnlist, NmodesDEIMAU, NmodesDEIMbU,"U_matrix"));
            // DEIMmatriceP = autoPtr<DEIM_functionP>( new DEIM_functionP(pEqnlist, NmodesDEIMAP, NmodesDEIMbP,"P_matrix"));
 
            // DEIMmatriceU->fieldsA = autoPtr<volVectorField>(new volVectorField(DEIMmatriceU->generateSubmeshMatrix(2, mesh, _U())));
            // DEIMmatriceU->fieldsB = autoPtr<volVectorField>(new volVectorField(DEIMmatriceU->generateSubmeshVector(2, mesh, _U())));
            // std::cout << "#####################################" << std::endl; 
            // DEIMmatriceP->fieldsA = autoPtr<volScalarField>(new volScalarField(DEIMmatriceP->generateSubmeshMatrix(2, mesh, _p())));
            // DEIMmatriceP->fieldsB = autoPtr<volScalarField>(new volScalarField(DEIMmatriceP->generateSubmeshVector(2, mesh, _p())));
            std::cout << "#####################################" << std::endl; 
            // // U Eqn
            // pFieldsA_U = DEIMmatriceU->generateSubFieldMatrix(_p());
            // UFieldsA_U = DEIMmatriceU->generateSubFieldMatrix(_U());
            // phiFieldsA_U = DEIMmatriceU->generateSubFieldMatrix(_phi());
            // pFieldsB_U = DEIMmatriceU->generateSubFieldVector(_p());
            // UFieldsB_U = DEIMmatriceU->generateSubFieldVector(_U());
            // phiFieldsB_U = DEIMmatriceU->generateSubFieldVector(_phi());


            // // P Eqn
            // pFieldsA_p = DEIMmatriceP->generateSubFieldMatrix(_p());
            // UFieldsA_p = DEIMmatriceP->generateSubFieldMatrix(_U());
            // phiFieldsA_p = DEIMmatriceP->generateSubFieldMatrix(_phi());
            // pFieldsB_p = DEIMmatriceP->generateSubFieldVector(_p());
            // UFieldsB_p = DEIMmatriceP->generateSubFieldVector(_U());
            // phiFieldsB_p = DEIMmatriceP->generateSubFieldVector(_phi());

            ModesUEig = Foam2Eigen::PtrList2Eigen(Umodes);
            ModesUEig.conservativeResize(ModesUEig.rows(), NmodesU);

            ModesPEig = Foam2Eigen::PtrList2Eigen(Pmodes);
            ModesPEig.conservativeResize(ModesPEig.rows(), NmodesP);

            // reduced operators
            ReducedMatricesA_U.resize(NmodesDEIMAU);
            ReducedVectorsB_U.resize(NmodesDEIMbU);

            ReducedMatricesA_p.resize(NmodesDEIMAP);
            ReducedVectorsB_p.resize(NmodesDEIMbP);

            //auto volU = ITHACAutilities::getMassMatrixFV(_U());
            //auto volP = ITHACAutilities::getMassMatrixFV(_p());
            //ReducedVectorsB_U = ModesUEig.transpose() * volU.asDiagonal()*DEIMmatriceU->MatrixOnlineB;

            // U
            for (int i = 0; i < NmodesDEIMAU; i++)
            {
                ReducedMatricesA_U[i] = ModesUEig.transpose() * DEIMmatriceU->MatrixOnlineA[i] * ModesUEig;
            }
            
            //ReducedMatricesA_U = ModesUEig.transpose() * volU.asDiagonal()* DEIMmatriceU->MatrixOnlineA * ModesUEig
            for (int i = 0; i < NmodesDEIMbU; i++)
            {
               ReducedVectorsB_U[i] = ModesUEig.transpose() * DEIMmatriceU->MatrixOnlineB;
            }


            // P
            for (int i = 0; i < NmodesDEIMAP; i++)
            {
                ReducedMatricesA_p[i] = ModesPEig.transpose() * DEIMmatriceP->MatrixOnlineA[i] * ModesPEig;
            }
            //ReducedMatricesA_p = ModesPEig.transpose() * volP.asDiagonal()* DEIMmatriceP->MatrixOnlineA * ModesPEig
            for (int i = 0; i < NmodesDEIMbP; i++)
            {
               ReducedVectorsB_p[i] = ModesPEig.transpose() * DEIMmatriceP->MatrixOnlineB;
            }
        

            std::cout << "###############" << ModesUEig.rows() << "###############" << std::endl;
            std::cout << "###############" << ModesUEig.cols() << "###############" << std::endl;
            std::cout << "###############" << DEIMmatriceU->MatrixOnlineA.size() << "###############" << std::endl;
        };

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
                sDRBMS.clear();
                argList& args = _args();
                Time& runTime = _runTime();
                runTime.setTime(0, 1);
                //meshPtr().movePoints(point0);    
                //pointField& pointOld = const_cast<pointField&> (meshPtr().oldPoints());
                //pointOld = point0;
                // _pimple.clear();
                Foam::dynamicFvMesh& mesh = meshPtr();
                _pimple = autoPtr<pimpleControl>
                               (
                                   new pimpleControl
                                   (
                                       mesh
                                   )
                           );

                
#include "createFields.H" 
        }
        // List<Eigen::MatrixXd> collocation(int Nmodes){
        //     List<Eigen::MatrixXd> collSys(2);
        //     List<label> magicPoints(DEIMmatrice->magicPoints());
        //     Eigen::MatrixXd A(magicPoints.size(), Nmodes);
        //     Eigen::MatrixXd b(magicPoints.size(), 1);

        //     for (label i = 0; i < Nmodes; i++)
        //     {
        //         volScalarField laplI(fvc::laplacian(nu(), Tmodes[i]));

        //         for (label j = 0; j < magicPoints.size(); j++)
        //         {
        //             A(j, i) =   laplI[magicPoints[j]];
        //             b(j)    =   S[magicPoints[j]];
        //         }
        //     }

        //     collSys[0] = A;
        //     collSys[1] = b;
        //     return collSys;
        // }

        };
        class ReducedFSIDEIM : public reductionProblem
        {
        public:
                explicit ReducedFSIDEIM(FSIDEIM& FoamProblem) : problem(&FoamProblem)
                {


                        for (int i = 0; i < problem->Umodes.size(); i++)
                        {
                            Umodes.append((problem->Umodes.toPtrList()[i]).clone());
                        }

                        for (int i = 0; i < problem->Pmodes.size(); i++)
                        {
                            Pmodes.append((problem->Pmodes.toPtrList()[i]).clone());
                        }

                }


                // Pointer to the full problem
                FSIDEIM* problem;
                // Online fields
                PtrList<volVectorField> Uonline;  
                PtrList<volScalarField> Ponline;
                PtrList<surfaceScalarField> Phionline;
                volScalarModes Pmodes;
                volVectorModes Umodes;
                scalar startTime = 0.0;
                scalar finalTime = 0.0;
                scalar timeStep = 0.0;
                scalar writeEvery = timeStep;
                scalar nextWrite = 0.0;
                void solveOnline(int NmodesUproj, int NmodesPproj, fileName folder = "./ITHACAoutput/Online/")
                {

                    Time& runTime = problem->_runTime();
                    surfaceScalarField& phi = problem->_phi();
                    dynamicFvMesh& mesh = problem->meshPtr();
                    fv::options& fvOptions = problem->_fvOptions();
                    pimpleControl& pimple = problem->_pimple();
                    volScalarField& p = problem->_p();
                    volVectorField& U = problem->_U();
                     /// constructing face velocity
                    //surfaceVectorField& Uf = *_Uf;
                    IOMRFZoneList& MRF = problem->_MRF();
                    singlePhaseTransportModel& laminarTransport = problem->_laminarTransport();
                    autoPtr<incompressible::turbulenceModel> turbulence = problem->turbulence;
                    //DEIM_functionU& DEIMmatriceU = problem->DEIMmatriceU();
                    //DEIM_functionP& DEIMmatriceP = problem->DEIMmatriceP();
                    label pRefCell = 0;
                    scalar pRefValue = 0.0;
                    instantList Times = runTime.times();
                    runTime.setEndTime(finalTime);
                    runTime.setTime(Times[1], 1);
                    runTime.setDeltaT(timeStep);
                    scalar nextWrite = startTime; 
                    label counter = 1;
                    turbulence->validate();

                    /// Exporting initial solution
                    ITHACAstream::exportSolution(U, name(counter), folder);
                    ITHACAstream::exportSolution(p, name(counter ), folder);
                    //ITHACAstream::exportSolution(phi, name(runTime.timeIndex()), folder);
                    ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
                    std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
                    Uonline.append(U.clone());
                    Ponline.append(p.clone());
                    //Phifield.append(phi.clone());
                    counter++;
                    nextWrite += writeEvery;

                    while (runTime.run())
                    {

                        runTime.setEndTime(finalTime);
                        runTime++;
                        Info<< "Time = " << runTime.timeName() << nl << endl;
                        p.storePrevIter();

                        MRF.correctBoundaryVelocity(U);
                        tmp<fvVectorMatrix> tUEqn
                        (
                            fvm::ddt(U) + fvm::div(phi, U)
                            + turbulence->divDevReff(U)
                        );
                        fvVectorMatrix& Ueqn1 = tUEqn.ref();
                        Ueqn1.relax();
                        fvOptions.constrain(Ueqn1);

                        volVectorField H = Ueqn1.H();
                        volScalarField A = Ueqn1.A();
                        List<Eigen::MatrixXd> RedLinSysU = Umodes.project(Ueqn1, NmodesUproj);
                        std::cout << "###############"  << "line 710" << "###############" << std::endl;
                        Eigen::MatrixXd thetaonAU = problem->DEIMmatriceU->onlineCoeffsAU(Ueqn1);
                        Eigen::MatrixXd thetaonBU = problem->DEIMmatriceU->onlineCoeffsBU(Ueqn1);
                        std::cout << "###############"  << "#####################" << "###############" << std::endl;
                        RedLinSysU[0] = EigenFunctions::MVproduct(problem->ReducedMatricesA_U, thetaonAU);

                        volVectorField gradpfull = -fvc::grad(p);
                        Eigen::MatrixXd projGrad = Umodes.project(gradpfull, NmodesUproj);
                        RedLinSysU[1] = EigenFunctions::MVproduct(problem->ReducedVectorsB_U, thetaonBU);
                        RedLinSysU[1] = RedLinSysU[1] + projGrad;

                        // a = reducedProblem::solveLinearSys(RedLinSysU, a, uresidual);
                        // Umodes.reconstruct(U, a, "U");
                        std::cout << "###############"  << "#####################" << "###############" << std::endl;

                        Eigen::VectorXd a = RedLinSysU[0].fullPivLu().solve(RedLinSysU[1]);
                        Eigen::VectorXd fullU = problem->ModesUEig* a;
                        volVectorField Ured("Ured", U);
                        U = Foam2Eigen::Eigen2field(Ured, fullU);
                       
                     ///////////////////////////////////////////////////////////////////////////////
                        volVectorField HbyA(constrainHbyA(1.0/A * H, U, p));
                        surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
                        List<Eigen::MatrixXd> RedLinSysP;
                       
                        while (pimple.correctNonOrthogonal())
                       {
                            
                            fvScalarMatrix peqn
                            (
                                fvm::laplacian(1/A, p) == fvc::div(phiHbyA)
                            );
                            peqn.setReference(pRefCell, pRefValue);

                            Eigen::MatrixXd thetaonAP = problem->DEIMmatriceP->onlineCoeffsAP(peqn);
                            Eigen::MatrixXd thetaonBP = problem->DEIMmatriceP->onlineCoeffsBP(peqn);
                            RedLinSysP = Pmodes.project(peqn, NmodesPproj);//???????
                            RedLinSysP[0] = EigenFunctions::MVproduct(problem->ReducedMatricesA_p, thetaonAP);
                            RedLinSysP[1]= EigenFunctions::MVproduct(problem->ReducedVectorsB_p, thetaonBP);
                            std::cout << "###############"  << "line 775" << "###############" << std::endl;
                            
                            //RedLinSysP[1] =  RedLinSysP[1].transpose();
                            //RedLinSysP[1]  = Eigen::Map<const Eigen::VectorXd>(RedLinSysP[1].data(), RedLinSysP[1].size());
                            //b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
                            //Pmodes.reconstruct(p, b, "p");
                            Eigen::MatrixXd b = RedLinSysP[0].fullPivLu().solve(RedLinSysP[1]);
                            Eigen::MatrixXd fullP = problem->ModesPEig* b; //reconstruct the solution
                            Eigen::MatrixXd fullPt = fullP.transpose();
                            Eigen::VectorXd v  = Eigen::Map<const Eigen::VectorXd>(fullPt.data(), fullPt.size());
                            volScalarField pred("pred", p);
                            p = Foam2Eigen::Eigen2field(pred, v);

                            if (pimple.finalNonOrthogonalIter())
                            {
                                phi = phiHbyA - peqn.flux();
                            }
                        }

                        p.relax();
                        U = HbyA - 1.0/A * fvc::grad(p);
                        U.correctBoundaryConditions();
                        if(problem->checkWrite(runTime))
                            {
                                ITHACAstream::exportSolution(U, name(counter), folder);
                                ITHACAstream::exportSolution(p, name(counter), folder);
                                Uonline.append(U.clone());
                                Ponline.append(p.clone());
                                //ITHACAstream::exportSolution(phi, name(runTime.timeIndex()), folder);
                                ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
                                std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
                                counter++;
                                nextWrite += writeEvery;
                            }   
                    }

                }// end solveonline
            };    


int main(int argc, char* argv[])
{
    
    // Construct the tutorial object
    FSIDEIM fsideim(argc, argv);
    FSIDEIM example(argc, argv);
    ITHACAparameters* para = ITHACAparameters::getInstance(fsideim.meshPtr(), fsideim._runTime());
    fsideim.startTime = 0;
    fsideim.finalTime = 2;
    fsideim.timeStep = 0.01; //0.01;
    fsideim.writeEvery = 0.1;
    //Perform the offline solve
    fsideim.offlineSolve();

    ITHACAPOD::getModes(fsideim.Ufield,   fsideim.Umodes,   fsideim._U().name(),   0, 0, 0, 10);
    ITHACAPOD::getModes(fsideim.Pfield,   fsideim.Pmodes,   fsideim._p().name(),   0, 0, 0, 10);
    //ITHACAPOD::getModes(fsideim.Ufield,   example.Umodes,   fsideim._U().name(),   0, 0, 0, 5);
    //ITHACAPOD::getModes(fsideim.Pfield,   example.Pmodes,   fsideim._p().name(),   0, 0, 0, 5);
    //ITHACAPOD::getModes(fsideim.Phifield, fsideim.Phimodes, fsideim._phi().name(), 0, 0, 0, 10);
  
    // Compute the offline part of the DEIM procedure
    fsideim.PODDEIM();
    fsideim.restart();
    //return 0;
    // ReducedFSIDEIM test(example);
    // test.startTime = 0;
    // test.finalTime = 2;
    // test.timeStep = 0.01; //0.01;
    // test.writeEvery = 0.1;


    Time& runTime = fsideim._runTime();
    surfaceScalarField& phi = fsideim._phi();
    dynamicFvMesh& mesh = fsideim.meshPtr();
    fv::options& fvOptions = fsideim._fvOptions();
    pimpleControl& pimple = fsideim._pimple();
    volScalarField& p = fsideim._p();
    volVectorField& U = fsideim._U();
     /// constructing face velocity
    //surfaceVectorField& Uf = *_Uf;
    IOMRFZoneList& MRF = fsideim._MRF();
    singlePhaseTransportModel& laminarTransport = fsideim._laminarTransport();
    autoPtr<incompressible::turbulenceModel> turbulence = fsideim.turbulence;
   
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    instantList Times = runTime.times();
    runTime.setEndTime(fsideim.finalTime);
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(fsideim.timeStep);
    scalar nextWrite = fsideim.startTime; 
    label counter = 1;
    
    turbulence->validate();
    fileName folder = "./ITHACAoutput/Online/";
    ITHACAstream::exportSolution(U, name(counter), folder);
    ITHACAstream::exportSolution(p, name(counter ), folder);
    //ITHACAstream::exportSolution(phi, name(runTime.timeIndex()), folder);
    ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
    std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
    //Ufield.append(U.clone());
    //Pfield.append(p.clone());
    counter++;
    nextWrite += fsideim.writeEvery;
    while (runTime.run())
    {

        runTime.setEndTime(fsideim.finalTime);
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;
        //p.storePrevIter();
        while (pimple.loop())
        {
                    MRF.correctBoundaryVelocity(U);
                    fvVectorMatrix Ueqn1
                    (
                        fvm::ddt(U) 
                        + fvm::div(phi, U) 
                        + MRF.DDt(U) 
                        + turbulence->divDevReff(U)
                        //== -fvc::grad(p)
                    );
                    //fvVectorMatrix& Ueqn1 = tUEqn.ref();
                    Ueqn1.relax();
                    fvOptions.constrain(Ueqn1);

                    Eigen::MatrixXd thetaonAU = fsideim.DEIMmatriceU->onlineCoeffsAU(Ueqn1);
                    Eigen::MatrixXd thetaonBU = fsideim.DEIMmatriceU->onlineCoeffsBU(Ueqn1);
                    //List<Eigen::MatrixXd> LinSysU;
                    List<Eigen::MatrixXd> LinSysU = (fsideim.Umodes).project(Ueqn1, 10);
std::cout << "=================== "<< "line 884"  <<"===================" << std::endl;
                    LinSysU[0] = EigenFunctions::MVproduct(fsideim.ReducedMatricesA_U, thetaonAU);
std::cout << LinSysU[0] << std::endl;

                    volVectorField gradpfull = -fvc::grad(p);
                    Eigen::MatrixXd projGrad = (fsideim.Umodes).project(gradpfull, 10);
                    //LinSysU[1] =  fsideim.ReducedVectorsB_U[0] * thetaonBU; 
                     LinSysU[1] = EigenFunctions::MVproduct(fsideim.ReducedVectorsB_U, thetaonBU);
                    LinSysU[1] = LinSysU[1] + projGrad;
std::cout << "###############"  << "line 898" << "###############" << std::endl;
                   Eigen::VectorXd a = LinSysU[0].fullPivLu().solve(LinSysU[1]);
                   Eigen::VectorXd fullU = fsideim.ModesUEig* a;
                   volVectorField Ured("Ured", U);
                   U = Foam2Eigen::Eigen2field(Ured, fullU);
std::cout << "###############"  << "line 903" << "###############" << std::endl;
                 ///////////////////////////////////////////////////////////////////////////////
                while(pimple.correct())
                {  
                    volVectorField H = Ueqn1.H();
                    volScalarField A = Ueqn1.A();
                    volVectorField HbyA(constrainHbyA(1.0/A * H, U, p));
                    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

                    List<Eigen::MatrixXd> LinSysP;
                    while (pimple.correctNonOrthogonal())
                    {
                        
                        fvScalarMatrix peqn
                        (
                            fvm::laplacian(1/A, p) == fvc::div(phiHbyA)
                        );
                        peqn.setReference(pRefCell, pRefValue);

                        Eigen::MatrixXd thetaonAP = fsideim.onlineCoeffsAP(peqn);
                        Eigen::MatrixXd thetaonBP = fsideim.onlineCoeffsBP(peqn);
                        
                        LinSysP = (fsideim.Pmodes).project(peqn, 10);
                        LinSysP[0] = EigenFunctions::MVproduct(fsideim.ReducedMatricesA_p, thetaonAP);
                        LinSysP[1]=  EigenFunctions::MVproduct(fsideim.ReducedVectorsB_p,  thetaonBP);

                        Eigen::MatrixXd b = LinSysP[0].fullPivLu().solve(LinSysP[1]);
std::cout << "###############"  << "line 930" << "###############" << std::endl;
                        Eigen::MatrixXd fullP = fsideim.ModesPEig* b; //reconstruct the solution
                        Eigen::MatrixXd fullPt = fullP.transpose();
                        Eigen::VectorXd v  = Eigen::Map<const Eigen::VectorXd>(fullPt.data(), fullPt.size());
                        volScalarField pred("pred", p);
                        p = Foam2Eigen::Eigen2field(pred, v);
                        //p = Foam2Eigen::Eigen2field(pred, fullP);
std::cout << "###############"  << "line 936" << "###############" << std::endl;
                        if (pimple.finalNonOrthogonalIter())
                        {
                            phi = phiHbyA - peqn.flux();
                        }
                    }

                    p.relax();
                    U = HbyA - 1.0/A * fvc::grad(p);
                    U.correctBoundaryConditions();
                }    
        } 
        //if(fsideim.checkWrite(runTime))
        //{
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);
            //Uonline.append(U.clone());
            //Ponline.append(p.clone());
            //ITHACAstream::exportSolution(phi, name(runTime.timeIndex()), folder);
            ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
            std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
            counter++;
            fsideim.nextWrite += fsideim.writeEvery;
        //}

    }
//test.solveOnline(5, 5);
    return 0;
}



 // Eigen::MatrixXd onlineCoeffsA(surfaceScalarField& fi)
        // {
        //     //Eigen::VectorXd fi2eigen = Foam2Eigen::field2Eigen(fi); 
        //     Eigen::MatrixXd theta(DEIMMomemtum->magicPointsAcol().size(), 1);
        //     fvVectorMatrix Aof = evaluate_momemtum(DEIMMomemtum->fieldA(), fi);
        //     Eigen::SparseMatrix<double> Mr;
        //     Eigen::VectorXd br;
        //     Foam2Eigen::fvMatrix2Eigen(Aof, Mr, br);
 
        //     for (int i = 0; i < DEIMMomemtum->magicPointsAcol().size(); i++)
        //     {
        //         int ind_row = DEIMMomemtum->localMagicPointsArow[i] + DEIMMomemtum->xyz_Arow()[i] *
        //                       DEIMMomemtum->fieldA().size();
        //         int ind_col = DEIMMomemtum->localMagicPointsAcol[i] + DEIMMomemtum->xyz_Acol()[i] *
        //                       DEIMMomemtum->fieldA().size();
        //         theta(i) = Mr.coeffRef(ind_row, ind_col);
        //     }
 
        //     return theta;
        // }




      // void PODDEIM()
        // {

        //     std::cout << "size of pEqnlist is: " << pEqnlist.size()  << std::endl;
        //     std::cout << "size of UEqnlist is: " << UEqnlist.size()  << std::endl;
        //     DEIMPpe = new DEIM_ppe(pEqnlist, 10, 1, "PMatrix");
        //     dynamicFvMesh& mesh = meshPtr();
        //     // Differential Operator
        //     DEIMPpe->fieldA = autoPtr<volScalarField>(new volScalarField(
        //                                 DEIMPpe->generateSubmeshMatrix(2, mesh, _p())));
        //     DEIMPpe->fieldB = autoPtr<volScalarField>(new volScalarField(
        //                                  DEIMPpe->generateSubmeshVector(2, mesh, _p())));
        //     // Source Terms
        //     ModesPEig = Foam2Eigen::PtrList2Eigen(Pmodes);
        //     ModesPEig.conservativeResize(ModesPEig.rows(), 5);
        //     ReducedMatricesP.resize(10);
        //     ReducedVectorsP.resize(1);

        //     for (int i = 0; i < 10; i++)
        //     {
        //         ReducedMatricesP[i] = ModesPEig.transpose() * DEIMPpe->MatrixOnlineA[i] * ModesPEig;
        //     }

        //     for (int i = 0; i < 1; i++)
        //     {
        //         ReducedVectorsP[i] = ModesPEig.transpose() * DEIMPpe->MatrixOnlineB;
        //     }
        // };


//class DEIM_DivTerm : public DEIM<fvVectorMatrix>
// {
//     public:
        
//         using DEIM::DEIM;

//         fvVectorMatrix evaluate_div(volVectorField& U, surfaceScalarField& phi)
//         {
            
//          tmp<fvVectorMatrix> UEqn1 = fvm::div(phi, U);
//           return UEqn1.ref();
//         }
//         autoPtr<volVectorField> fieldA;
//         autoPtr<volVectorField> fieldB;

// };

// class DEIM_ppe : public DEIM<fvScalarMatrix>
// {
//     public:
//         //friend class FSIDEIM;
//         using DEIM::DEIM;

//         DEIM_ppe(PtrList<fvScalarMatrix>& s, label MaxModesA, label MaxModesB, word MatrixName)
//         : DEIM(s, MaxModesA, MaxModesB, MatrixName) 
//         {
//             std::cout << "########## DEIM for the pressure eqn ##############" << std::endl;
//         }

//         fvScalarMatrix evaluate_ppe(fvVectorMatrix& UEqn, volScalarField& p, surfaceScalarField& fi)
//         {
            
//             volScalarField rAU(1.0/ UEqn.A());
//             tmp<volScalarField> rAtU(rAU);

//             fvScalarMatrix pEqn( fvm::laplacian(rAtU(), p) == fvc::div(fi) );
//             return pEqn;

//         } 
//         autoPtr<volScalarField> fieldA;
//         autoPtr<volScalarField> fieldB;   
// };
// class DEIM_test : public DEIM<volScalarField>
// {
//     public:

//         using DEIM::DEIM;

//         autoPtr<volScalarField> fieldA;
//         autoPtr<volScalarField> fieldB;   
// };

// class DEIM_Flux : public DEIM<surfaceScalarField>
// {
//     public:

//         using DEIM::DEIM;
//         surfaceScalarField evaluate_flux(fvVectorMatrix& UEqn, volScalarField& p, volVectorField& U)
//         {
            
//             volScalarField rAU(1.0/UEqn.A());
//             volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
//             surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));


//             return phiHbyA;
//         } 
//         autoPtr<surfaceScalarField> fields;
//         autoPtr<surfaceScalarField> subfield; 
// };


        // fvVectorMatrix evaluate(volVectorField& U)
        // {
            
        //     fvVectorMatrix UEqn2(fvm::ddt(U)+ turbulence->divDevReff(U));
        //     return UEqn2;
        // };

        // fvVectorMatrix evaluate_div(volVectorField& U, surfaceScalarField& phi)
        // {
            
        //  fvVectorMatrix UEqn1(fvm::div(phi, U));
        //   return UEqn1;
        // }
        // surfaceScalarField evaluate_flux(fvVectorMatrix& Ueqn, volScalarField& p, volVectorField& U)
        // {
            
        //     volScalarField rAU(1.0/Ueqn.A());
        //     volVectorField HbyA(constrainHbyA(rAU*Ueqn.H(), U, p));
        //     surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));


        //     return phiHbyA;
        // } 
        
        // fvScalarMatrix evaluate_ppe(fvVectorMatrix& Ueqn, volScalarField& p, surfaceScalarField& fi)
        // {
            
        //     volScalarField rAU(1.0/ Ueqn.A());
        //     //tmp<volScalarField> rAtU(rAU);

        //     fvScalarMatrix peqn( fvm::laplacian(rAU, p) == fvc::div(fi) );
        //     return peqn;

        // }  

