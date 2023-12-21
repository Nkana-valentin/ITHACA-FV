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
 \*---------------------------------------------------------------------------*/
std::cout << "*****create fields of a  dynamic case****\n"<< std::endl; 
Info << "Reading field p\n" << endl;
_p = autoPtr<volScalarField>
     (
         new volScalarField
         (
             IOobject
             (
                 "p",
                 runTime.timeName(),
                 mesh,
                 IOobject::MUST_READ,
                 IOobject::AUTO_WRITE
             ),
             mesh
         )
     );
volScalarField& p = _p();

Info << "Reading field U\n" << endl;
_U = autoPtr<volVectorField>
     (
         new volVectorField
         (
             IOobject
             (
                 "U",
                 runTime.timeName(),
                 mesh,
                 IOobject::MUST_READ,
                 IOobject::AUTO_WRITE
             ),
             mesh
         )
     );
volVectorField& U = _U();

Info << "Reading/calculating face flux field phi\n" << endl;

_phi = autoPtr<surfaceScalarField>
       (
           new surfaceScalarField
           (
               IOobject
               (
                   "phi",
                   runTime.timeName(),
                   mesh,
                   IOobject::READ_IF_PRESENT,
                   IOobject::AUTO_WRITE
               ),
               linearInterpolate(U) & mesh.Sf()
           )
       );
       
surfaceScalarField& phi = _phi();


_Uf = autoPtr<surfaceVectorField>
       (
           new surfaceVectorField
           (
               IOobject
               (
                   "Uf",
                   runTime.timeName(),
                   mesh,
                   IOobject::READ_IF_PRESENT,
                   IOobject::AUTO_WRITE
               ),
               fvc::interpolate(U)
           )
       );
       
surfaceVectorField& Uf = _Uf();

_T = autoPtr<volVectorField>
     (
         new volVectorField
         (
             IOobject
             (
                 "T",
                 runTime.timeName(),
                 mesh,
                 IOobject::MUST_READ,
                 IOobject::AUTO_WRITE
             ),
             mesh
         )
     );
volVectorField& T = _T();

const pointMesh& pMesh = pointMesh::New(mesh);
_pd = autoPtr<pointVectorField>
       (
           new pointVectorField
           (
               IOobject
                (
                    "pointDisplacement",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                pMesh,
                dimensionedVector(dimLength, Zero),
                pointPatchVectorField::calculatedType()
           )
       );
pointVectorField& pd = _pd();

label pRefCell = 0;
scalar pRefValue = 0.0;

setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
_laminarTransport = autoPtr<singlePhaseTransportModel>
                    (
                        new singlePhaseTransportModel( U, phi )
                    );
singlePhaseTransportModel& laminarTransport = _laminarTransport();
//std::cout << "##############createFields line 114 #################\n"<< std::endl; 

turbulence = autoPtr<incompressible::turbulenceModel>
             (
                 incompressible::turbulenceModel::New(U, phi, laminarTransport)
             );
//#include "createMRF.H"
_MRF = autoPtr<IOMRFZoneList>
       (
           new IOMRFZoneList(mesh)
       );

_fvOptions = autoPtr<fv::options>(new fv::options(mesh));


    dyndict = new IOdictionary
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
    sDRBMS = autoPtr<sixDoFRigidBodyMotionSolver>
    (
    new sixDoFRigidBodyMotionSolver(mesh, *dyndict)
    );
//std::cout << "======================= line 162    ===================================" << std::endl;
    pimpleControl& pimple = _pimple();
    //std::cout << "======================= line 164    ===================================" << std::endl;
    correctPhi = pimple.dict().getOrDefault("correctPhi",mesh.dynamic());
    //std::cout << "======================= line 166    ===================================" << std::endl;
    checkMeshCourantNo = pimple.dict().getOrDefault("checkMeshCourantNo",false);
    //std::cout << "======================= line 167    ===================================" << std::endl;
    moveMeshOuterCorrectors = pimple.dict().getOrDefault("moveMeshOuterCorrectors",false);
    //std::cout << "======================= line 170    ===================================" << std::endl;
