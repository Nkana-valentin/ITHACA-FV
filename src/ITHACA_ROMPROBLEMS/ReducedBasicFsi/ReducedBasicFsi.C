#include "ReducedBasicFsi.H"


/// Constructor Null
ReducedBasicFsi::ReducedBasicFsi() {}

ReducedBasicFsi::ReducedBasicFsi(fsiBasic& FoamPb): problem(&FoamPb)
        
{
            for (int i = 0; i < problem->inletIndex.rows(); i++)
            {
                Umodes.append((problem->liftfield[i]).clone());
            }
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
            
  
}

void ReducedBasicFsi::solveOnline_Pimple(scalar mu_now, int NmodesUproj, int NmodesPproj, int NmodesDproj, fileName folder)
    {
        

        Eigen::VectorXd uresidual = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::VectorXd presidual = Eigen::VectorXd::Zero(NmodesPproj);

        Eigen::MatrixXd a = Eigen::VectorXd::Zero(NmodesUproj);
        Eigen::MatrixXd b = Eigen::VectorXd::Zero(NmodesPproj);
         Eigen::MatrixXd bOld = b;
        //Eigen::MatrixXd c = Eigen::VectorXd::Zero(NmodesDproj);
        Time& runTime = problem->_runTime();
        surfaceScalarField& phi = problem->_phi();
        dynamicFvMesh& mesh = problem->meshPtr();
        const pointMesh& pMesh = pointMesh::New(mesh);
        fv::options& fvOptions = problem->_fvOptions();
        pimpleControl& pimple = problem->_pimple();
        volScalarField& p = problem->_p();
        volVectorField& U = problem->_U();
        pointVectorField& pointDisplacement= problem->_pointDisplacement();
        // ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
        // ITHACAstream::read_fields(Pfield, p, "./ITHACAoutput/Offline/");
        // U = Ufield[0];
        // p = Pfield[0];
        // pointVectorField& pointDisplacement 
        // = const_cast<pointVectorField&>
        // (mesh.lookupObject<pointVectorField>("pointDisplacement"));
        //ITHACAstream::read_fields(Dfield, pointDisplacement, "./ITHACAoutput/Offline/");
        //sixDoFRigidBodyMotionSolver& sDRBMS = problem->sDRBMS();
        //To set initial condition properly
        // pointDisplacement.primitiveFieldRef() = problem->Dfield[0].primitiveFieldRef(); 
        // problem->sDRBMS().pointDisplacement().primitiveFieldRef() = problem->Dfield[0].primitiveField();
        // problem->sDRBMS().pointDisplacement().primitiveFieldRef() = pointDisplacement.primitiveFieldRef();

        IOMRFZoneList& MRF = problem->_MRF();
        singlePhaseTransportModel& laminarTransport = problem->_laminarTransport();
        autoPtr<incompressible::turbulenceModel> turbulence = problem->turbulence;
        instantList Times = runTime.times();
        runTime.setEndTime(finalTime);
        // Declare modal coefficients for velocity,  pressure and pointDisplacement
        //Initial conditions
        //a = ITHACAutilities::getCoeffs(U, Umodes, NmodesUproj, true);
        //b = ITHACAutilities::getCoeffs(p, Pmodes, NmodesPproj, true);//????
        //c = ITHACAutilities::getCoeffs(pointDisplacement, Dmodes, NmodesDproj, false);
////////////////////////////////////////////////////////////////
        //Umodes.reconstruct(U, a, "U");
       // Pmodes.reconstruct(p, b, "p");
        //Dmodes.reconstruct(pointDisplacement, c, "pointDisplacement");
// To solve the rbfi system to obtain the weights
PodIpointDispl(problem->coeffL2, problem->CylDispl, NmodesDproj);
        
    
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

        Eigen::MatrixXd pdCoeff;
        pdCoeff.resize(NmodesDproj, 1);
        //- Current time index (used for updating)
        label curTimeIndex_ = -1;
        pointField points0 = mesh.points();
        Eigen::MatrixXd muEval;
        muEval.resize(1, 1);


        // PIMPLE algorithm starts here
        Info<< "\nStarting time loop\n" << endl;
        while (runTime.run())
        {

            runTime.setEndTime(finalTime);
            runTime++;
            //p.storePrevIter();
            Info << "Time = " << runTime.timeName() << nl << endl;

            while (pimple.loop())
            {
                if (pimple.firstIter() || moveMeshOuterCorrectors)
                {
                  
                    //- Store the motion state at the beginning of the time-step
                    const Time& t = mesh.time();
                    bool firstIter = false;
                    if (curTimeIndex_ != mesh.time().timeIndex())
                    {
                        sDRBM.newTime();
                        curTimeIndex_ = mesh.time().timeIndex();
                        firstIter = true;
                    }
                    Foam::dimensionedVector g("g", dimAcceleration, Zero);
                    if (mesh.time().foundObject<uniformDimensionedVectorField>("g"))
                    {
                        g = mesh.time().lookupObject<uniformDimensionedVectorField>("g");
                    }
                    else
                    {
                        dictCoeffs.readIfPresent("g", g);
                    }
                    const scalar ramp = 1.0;
                    dictionary forcesDict;
                    forcesDict.add("type", functionObjects::forces::typeName);
                    forcesDict.add("patches", dictCoeffs.get<wordRes>("patches"));
                    forcesDict.add("rhoInf", 1.0);
                    forcesDict.add("rho", dictCoeffs.getOrDefault<word>("rho", "rho"));
                    forcesDict.add("CofR",sDRBM.centreOfRotation());
                    
                    Foam::functionObjects::forces romforces("romforces", mesh, forcesDict);
                    romforces.calcForcesMoment(); //calcForcesMoments()
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
                         t.deltaTValue(),
                         t.deltaT0Value()
                    );
                    // To append the linear velocities
                    //velx.append(sDRBM.v().x());
                    vely.append(sDRBM.v().y());
                    //for (int i = 0; i < NmodesDproj; i++)
                    //{
                        //New value of the parameter
                        muEval(0, 0)=  sDRBM.centreOfMass().y(); 
                        //pdCoeff(i, 0) = problem->rbfSplines[i]->eval(muEval);
                        pdCoeff(0, 0) = problem->rbfSplines[0]->eval(muEval);
std::cout << "=================" << pdCoeff(0,0) << "===============" << std::endl;
                        // RedPdCoeff.append(pdCoeff(0, 0));
                    //}
                    // Reconstruction of the pointdisplacement
                    Dmodes.reconstruct(pointDisplacement, pdCoeff, "pointDisplacement");
                    // update the new pointDisplacement                 
                    problem->sDRBMS().pointDisplacement().primitiveFieldRef() 
                    = pointDisplacement.primitiveFieldRef();                  
                    // Displacement has changed. Update boundary conditions // 
                    pointConstraints::New
                    (
                       pointDisplacement.mesh()

                    ).constrainDisplacement(problem->sDRBMS().pointDisplacement()); //????????????

                    mesh.movePoints(problem->sDRBMS().curPoints());

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
                    fvm::ddt(U) 
                    + fvm::div(phi, U)
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
                    bOld = b;
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
                        //b = bOld + mesh.fieldRelaxationFactor("p") * (b - bOld);
                        //b = bOld + b;
                        Pmodes.reconstruct(p, b, "p");

                        if (pimple.finalNonOrthogonalIter())
                        {
                            phi = phiHbyA - pEqn.flux();
                        }
                    }
                    // Explicitly relax pressure for momentum corrector
                    p.relax();
                    //b = bOld + mesh.fieldRelaxationFactor("p") * (b - bOld);
                    //Pmodes.reconstruct(p, b, "p");
                    U = HbyA - rAtU * fvc::grad(p); //p
                    U.correctBoundaryConditions();
                    fvOptions.correct(U);
                    // Correct Uf if the mesh is moving
                    fvc::correctUf(Uf, U, phi);

                    // Make the fluxes relative to the mesh motion
                    fvc::makeRelative(phi, U);

                }// end of the pimple.correct()
                //Umodes.reconstruct(U, a, "U"); 
                //Pmodes.reconstruct(p, b, "p");

            }// end of the pimple.loop()
            if(checkWrite(runTime))
            {
                // romforcex.append(romforces.forceEff().x());
                // romforcey.append(romforces.forceEff().y()); 

                // centerofmassx.append(sDRBM.centreOfMass().x());
                centerofmassy.append(sDRBM.centreOfMass().y());
                // centerofmassz.append(sDRBM.centreOfMass().z());
                // // To append the linear velocities
                // velx.append(sDRBM.v().x());
                // vely.append(sDRBM.v().y());
                // velz.append(sDRBM.v().z());
                pdcoeffrbf.append(pdCoeff(0, 0));
                CoeffP.append(b);
                CoeffU.append(a);

                ITHACAstream::exportSolution(U, name(counter), folder);
                ITHACAstream::exportSolution(p, name(counter), folder);
                ITHACAstream::exportSolution(pointDisplacement, name(counter), folder);
                ITHACAstream::writePoints(mesh.points(), folder, name(counter) + "/polyMesh/");
                std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
                UredFields.append(U.clone());
                PredFields.append(p.clone());
                DfieldRbf.append(pointDisplacement.clone());
                counter++;
                nextWrite += writeEvery;
            }

        }

    }

    bool ReducedBasicFsi::checkWrite(Time& timeObject)
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


    void ReducedBasicFsi::PodIpointDispl(Eigen::MatrixXd coeffL2, Eigen::MatrixXd muu,  label NPdModes)
    {
            if (NPdModes == 0)
            {
                NPdModes = Dmodes.size();
            }

            problem->samples.resize(NPdModes);
            problem->rbfSplines.resize(NPdModes);
            Eigen::MatrixXd weights;
          
            for (label i = 0; i < NPdModes; i++) // i is the nnumber of th mode
            {
                word weightName = "wRBF_M" + name(i + 1);

                if (ITHACAutilities::check_file("./ITHACAoutput/weights/" + weightName))
                {
                    problem->samples[i] = new SPLINTER::DataTable(1, 1);
                    for (label j = 0; j < coeffL2.cols(); j++) // j is the number of the nut snapshot
                    {
                        problem->samples[i]->addSample(muu.row(j), coeffL2(i, j));
                    }

                    ITHACAstream::ReadDenseMatrix(weights, "./ITHACAoutput/weights/", weightName);
                    problem->rbfSplines[i] = new SPLINTER::RBFSpline(*(problem->samples)[i],
                                                            SPLINTER::RadialBasisFunctionType::GAUSSIAN, weights);
                    std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
                } 
               
                else
                {
                    problem->samples[i] = new SPLINTER::DataTable(1, 1);

                    for (label j = 0; j < coeffL2.cols();j++) // j is the number of the nut snapshot
                    {
                        problem->samples[i]->addSample(muu.row(j), coeffL2(i, j));
                    }
                    problem->rbfSplines[i] = new SPLINTER::RBFSpline(*(problem->samples)[i],
                                                            SPLINTER::RadialBasisFunctionType::GAUSSIAN);
                    ITHACAstream::SaveDenseMatrix(problem->rbfSplines[i]->weights,
                                                  "./ITHACAoutput/weights/", weightName);
                    std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
                }
            }
        }