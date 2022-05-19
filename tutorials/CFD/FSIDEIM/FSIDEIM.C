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
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "sixDoFRigidBodyMotionSolver.H"
#include "ITHACAutilities.H"
#include "Foam2Eigen.H"
#include "DEIM.H"
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "IOmanip.H"
#include <chrono>
#include "Modes.H"
#include "reductionProblem.H"
#include "ReducedProblem.H"
#include "EigenFunctions.H"
#include <Eigen/SVD>
#include <Eigen/SparseLU>
#include "GeometricFields.H"
#include "transformGeometricField.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>


class DEIM_DivTerm : public DEIM<fvVectorMatrix>
{
    public:
        
        using DEIM::DEIM;

        fvVectorMatrix evaluate_div(volVectorField& U, surfaceScalarField& phi)
        {
            
         
          return fvm::div(phi, U);
        }
        autoPtr<volVectorField> fieldA;
        autoPtr<volVectorField> fieldB;

};

class DEIM_ppe : public DEIM<fvScalarMatrix>
{
    public:
        //friend class FSIDEIM;
        using DEIM::DEIM;

        DEIM_ppe(PtrList<fvScalarMatrix>& s, label MaxModesA, label MaxModesB, word MatrixName)
        : DEIM(s, MaxModesA, MaxModesB, MatrixName) 
        {
            std::cout << "########## DEIM for the pressure eqn ##############" << std::endl;
        }

        static fvScalarMatrix evaluate_ppe(fvVectorMatrix& UEqn, volScalarField& p, surfaceScalarField& fi)
        {
            
            volScalarField rAU(1.0/ UEqn.A());
            tmp<volScalarField> rAtU(rAU);

            fvScalarMatrix pEqn( fvm::laplacian(rAtU(), p) == fvc::div(fi) );
            return pEqn;

        } 
        autoPtr<volScalarField> fieldA;
        autoPtr<volScalarField> fieldB;   
};
class DEIM_test : public DEIM<volScalarField>
{
    public:

        using DEIM::DEIM;

        autoPtr<volScalarField> fieldA;
        autoPtr<volScalarField> fieldB;   
};

class DEIM_Flux : public DEIM<surfaceScalarField>
{
    public:

        using DEIM::DEIM;
        static surfaceScalarField evaluate_flux(fvVectorMatrix& UEqn, volScalarField& p, volVectorField& U)
        {
            
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));


            return phiHbyA;
        } 
        autoPtr<surfaceScalarField> fields;
        autoPtr<surfaceScalarField> subfield; 
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
        scalar writeEvery = 0.0;
        scalar nextWrite = 0.0;
        autoPtr<argList> _args;
        autoPtr<pimpleControl> _pimple;
        autoPtr<fv::options> _fvOptions;
        autoPtr<Time> _runTime;

        PtrList<volScalarField> Pfield;
        PtrList<volVectorField> Ufield;
        PtrList<surfaceScalarField> Phifield;
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
      
        autoPtr<incompressible::turbulenceModel> turbulence;
        autoPtr<singlePhaseTransportModel> _laminarTransport;
        autoPtr<sixDoFRigidBodyMotionSolver> sDRBMS;
        IOdictionary* dyndict;
      
        /// DEIM fields
        DEIM_Flux* FluxObject = nullptr;
        DEIM_DivTerm* DEIMObject = nullptr;
        DEIM_ppe* DEIMPpe =nullptr;
        //autoPtr<DEIM_ppe> DEIMPpe;
        PtrList<fvScalarMatrix> pEqnlist;
        PtrList<fvVectorMatrix> UEqnlist;
        PtrList<fvVectorMatrix> DivList;
        PtrList<surfaceScalarField> Philist;
        // std::vector<Eigen::MatrixXd> ReducedMatricesU;
        // std::vector<Eigen::MatrixXd> ReducedVectorsU;
        // std::vector<Eigen::MatrixXd> ReducedMatricesP;
        // std::vector<Eigen::MatrixXd> ReducedVectorsP;
        // Eigen::MatrixXd ModesUEig;
        // Eigen::MatrixXd ModesPEig;
        // Eigen::MatrixXd ModesPhiEig;
        int NmodesU;
        int NmodesDEIMA; 
        int NmodesDEIMB;
        //autoPtr<incompressible::turbulenceModel> turbulence;
        // List of Modes for P
        volScalarModes Pmodes;
        // List of Modes for U
        volVectorModes Umodes;
        // List of Modes for Phi
        surfaceScalarModes Phimodes;
        
        void offlineSolve(word folder = "./ITHACAoutput/Offline/") 
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
            //IOMRFZoneList& MRF = _MRF();
            singlePhaseTransportModel& laminarTransport = _laminarTransport();
            turbulence = autoPtr<incompressible::turbulenceModel>
            (
                     incompressible::turbulenceModel::New(U, phi, laminarTransport)
            );
            label pRefCell = 0;
            scalar pRefValue = 0.0;
            instantList Times = runTime.times();
            runTime.setEndTime(finalTime);
            runTime.setTime(Times[1], 1);
            runTime.setDeltaT(timeStep);
            nextWrite = startTime; 
        
            ITHACAstream::exportSolution(U, name(runTime.timeIndex()), folder);
            ITHACAstream::exportSolution(p, name(runTime.timeIndex()), folder);
            ITHACAstream::exportSolution(phi, name(runTime.timeIndex()), folder);
            ITHACAstream::writePoints(mesh.points(), folder, name(runTime.timeIndex()) + "/polyMesh/");
            Ufield.append(U.clone());
            Pfield.append(p.clone());
            Phifield.append(phi.clone());
  
            while(runTime.run()){

                runTime.setEndTime(finalTime);
                ++runTime;
                Info<< "Time = " << runTime.timeName() << nl << endl;
                while(pimple.loop())
                {
                    ///From UEqn.H
                    fvVectorMatrix divterm = DEIMObject->evaluate_div(U, phi);
                    //fvVectorMatrix UEqn = fvm::ddt(U) + DEIMObject->evaluate_div(U, phi) + turbulence->divDevReff(U);
                    fvVectorMatrix UEqn = fvm::ddt(U) + divterm + turbulence->divDevReff(U);
                    UEqn.relax();
                    fvOptions.constrain(UEqn);
                    //if (pimple.momentumPredictor())
                   // {
                       solve(UEqn == -fvc::grad(p));
                       fvOptions.correct(U);
                    //}   

                    /// From pEqn.H
                    volScalarField A_inv(1.0/ UEqn.A());
                    volVectorField H = UEqn.H();
                    //volVectorField HbyA(constrainHbyA(A_inv*UEqn.H(), U, p));
                    // Computing HbyA field = H/A for ease of calculation
                    volVectorField HbyA = A_inv * H;
                    UEqnlist.append(UEqn.clone());
                    DivList.append(divterm.clone());
                    // Evaluate the flux
                    surfaceScalarField Flux  = FluxObject->evaluate_flux(UEqn, p, U);
                    Philist.append(Flux.clone()); //OK
                    // Forming the pressure correction equation:
                    // Nab(A^-1 Nab(p)) = Nab.(A^-1 * H)
                    // The LHS can be defined using the laplacian operator in OpenFOAM as:
                    while (pimple.correctNonOrthogonal())
                    {
                        fvScalarMatrix pEqn = DEIMPpe->evaluate_ppe(UEqn, p, Flux);
                        // Solving the pressure correction equation
                         pEqnlist.append(pEqn.clone());
                        //Info << pEqn << endl;
                        pEqn.setReference(pRefCell, pRefValue);
                         // Solving the pressure correction equation.
                        pEqn.solve();
                     
                       
                    }    
                    // Under-relaxing the pressure equation using explicit relaxation:
                    //p = alpha*p + (1.0 - alpha)*p_old;
                    /// Update the new flux
                    //phi = phiHbyA - pEqn.flux();
                    p.relax();
                    U = HbyA - A_inv*fvc::grad(p);
                    // Updating the flux field with newly updated velocity field.
                    phi = linearInterpolate(U) & mesh.Sf();//???????
                    // Updating boundary conditions for both p and U fields.
                    U.correctBoundaryConditions();
                    fvOptions.correct(U);
                    //p.correctBoundaryConditions();
                    // Updating old_pressure field with new values
                    //p_old = p;

                    // Correct Uf if the mesh is moving
                    //fvc::correctUf(_Uf, U, phi);
                     
                    // Make the fluxes relative to the mesh motion
                    //fvc::makeRelative(phi, U);
                }        
              
                ITHACAstream::exportSolution(U, name(runTime.timeIndex()), folder);
                ITHACAstream::exportSolution(p, name(runTime.timeIndex()), folder);
                ITHACAstream::exportSolution(phi, name(runTime.timeIndex()), folder);
                ITHACAstream::writePoints(mesh.points(), folder, name(runTime.timeIndex()) + "/polyMesh/");
                Ufield.append(U.clone());
                Pfield.append(p.clone());
                Phifield.append(phi.clone());
               
               
            }

        };


        void PODDEIM()
        {

            //std::cout << "size of pEqnlist is: " << pEqnlist.size()  << std::endl;
            DEIMPpe = new DEIM_ppe(pEqnlist,2,2, "PMatrix");

        }
        

        fvVectorMatrix evaluate(volVectorField& U)
        {
            
            fvVectorMatrix UEqn2( fvm::ddt(U)+ turbulence->divDevReff(U));
            return UEqn2;
        };
    
};


int main(int argc, char* argv[])
{
    

    FSIDEIM fsideim(argc, argv);
    ITHACAparameters* para = ITHACAparameters::getInstance(fsideim.meshPtr(), fsideim._runTime());
    fsideim.startTime = 0;
    fsideim.finalTime = 1;
    fsideim.timeStep = 0.01; //0.01;
    fsideim.writeEvery = 0.01;
    fsideim.offlineSolve();

    ITHACAPOD::getModes(fsideim.Ufield, fsideim.Umodes, fsideim._U().name(), fsideim.podex, 0, 0, 10);
    ITHACAPOD::getModes(fsideim.Pfield, fsideim.Pmodes, fsideim._p().name(), fsideim.podex, 0, 0, 10);
    ITHACAPOD::getModes(fsideim.Phifield, fsideim.Phimodes, fsideim._phi().name(), fsideim.podex, 0, 0, 10);
    // fsideim.DEIMPpe->fieldA = autoPtr<volScalarField>(new volScalarField(
    //                                       fsideim.DEIMPpe->generateSubmeshMatrix(2, fsideim.meshPtr(), fsideim.p)));
    fsideim.PODDEIM();
    //DEIM_momemtum* DEIMDiv = new DEIM_momemtum(fsideim.UEqnlist, 2, 2, "UMatrix");
    //std::cerr << "gfgfghfghhhhhhhhhhhhhhhhhhhhhhhhhhhhhh" << std::endl;
    //DEIM_ppe* DEIMPpe = new DEIM_ppe(fsideim.pEqnlist, 10, 1, "PMatrix");
    //DEIM_test test(DEIM_test(fsideim.Pfield, 2, "funcname", fsideim._p().name()));
    return 0;
}




 // void PODDEIM()
 //        {

 //            std::cout << "size of pEqnlist is: " << pEqnlist.size()  << std::endl;
 //            DEIMPpe = new DEIM_ppe(pEqnlist, NmodesDEIMA, NmodesDEIMB, "PMatrix");
 //              fvMesh& mesh  =  const_cast<fvMesh&>(T.mesh());
 //            //     // fvMesh& mesh  =  const_cast<fvMesh&>(U.mesh());
 //            //     // // Differential Operator
 //            //     // DEIMmatrice_U->fieldA = autoPtr<volVectorField>(new volVectorField(
 //            //     //                             DEIMmatrice_M->generateSubmeshMatrix(2, mesh, U)));
 //            //     // DEIMmatrice_U->fieldB = autoPtr<volVectorField>(new volVectorField(
 //            //     //                             DEIMmatrice_M->generateSubmeshVector(2, mesh, U)));
 //            //         // // Source Terms
 //            //         // ModesUEig = Foam2Eigen::PtrList2Eigen(Tmodes);
 //            //         // ModesUEig.conservativeResize(ModesUEig.rows(), NmodesU);
 //            //         // ReducedMatricesA.resize(NmodesDEIMA);
 //            //         // ReducedVectorsB.resize(NmodesDEIMB);

 //            //         // for (int i = 0; i < NmodesDEIMA; i++)
 //            //         // {
 //            //         //     ReducedMatricesA[i] = ModesUEig.transpose() * DEIMmatrice_M->MatrixOnlineA[i] *
 //            //         //                           ModesUEig;
 //            //         // }

 //            //         // for (int i = 0; i < NmodesDEIMB; i++)
 //            //         // {
 //            //         //     ReducedVectorsU[i] = ModesUEig.transpose() * DEIMmatrice_M->MatrixOnlineB;
 //            //         // }
 //        };



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
