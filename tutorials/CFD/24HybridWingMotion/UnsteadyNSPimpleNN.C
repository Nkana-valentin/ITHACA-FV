#include "UnsteadyNSPimpleNN.H"

UnsteadyNSPimpleNN::UnsteadyNSPimpleNN() {}
UnsteadyNSPimpleNN::UnsteadyNSPimpleNN(int argc, char* argv[])
            :
        fsiBasic(argc, argv)
        {

            _nut = autoPtr<volScalarField>
               (
                   new volScalarField
                   (
                       IOobject
                       (
                           "nut",
                           _runTime().timeName(),
                           meshPtr(),
                           IOobject::MUST_READ,
                           IOobject::AUTO_WRITE
                       ),
                       meshPtr()
                   )
               );
        }

        UnsteadyNSPimpleNN::~UnsteadyNSPimpleNN(){}
     

    void UnsteadyNSPimpleNN::loadNet(word filename)
    {
        std::string Msg = filename +
                          " is not existing, please run the training stage of the net with the correct number of modes for U and Nut";
        M_Assert(ITHACAutilities::check_file(filename), Msg.c_str());
        netTorchscript = torch::jit::load(filename);
        cnpy::load(bias_inp, "ITHACAoutput/NN/minInp_" + name(NUmodes) + "_" + name(NNutModes) + ".npy");
        cnpy::load(scale_inp, "ITHACAoutput/NN/scaleInp_" + name(NUmodes) + "_" + name(NNutModes) + ".npy");
        cnpy::load(bias_out, "ITHACAoutput/NN/minOut_" + name(NUmodes) + "_" + name(NNutModes) + ".npy");
        cnpy::load(scale_out, "ITHACAoutput/NN/scaleOut_" + name(NUmodes) + "_" + name(NNutModes) + ".npy");
        netTorchscript.eval();
    }

    void UnsteadyNSPimpleNN::loadLstmNet(word filename)
    {
        std::string Msg = filename + " is not existing, please run the training stage of  the Lstm for Nut";
        M_Assert(ITHACAutilities::check_file(filename), Msg.c_str());
        lstmTorchscript = torch::jit::load(filename);
        cnpy::load(bias_inp, "ITHACAoutput/NN/LstmMinInp_"  +   name(NNutModes) + ".npy");
        cnpy::load(scale_inp, "ITHACAoutput/NN/LstmScaleInp_"   + name(NNutModes) + ".npy");
        cnpy::load(bias_out, "ITHACAoutput/NN/LstmMinOut_"     +  name(NNutModes) + ".npy");
        cnpy::load(scale_out, "ITHACAoutput/NN/LstmScaleOut_"    +  name(NNutModes) + ".npy");
        lstmTorchscript.eval();
    }

    // This function computes the coefficients which are later used for training
    void UnsteadyNSPimpleNN::getTurbNN()
    {
        if (!ITHACAutilities::check_folder("ITHACAoutput/NN/coeffs"))
        {
            mkDir("ITHACAoutput/NN/coeffs");
            /// Read Fields for Train
            PtrList<volVectorField> UfieldData;
            PtrList<volScalarField> PfieldData;
            PtrList<volScalarField> nutFieldsData;
            ITHACAstream::read_fields(UfieldData, _U(), "./ITHACAoutput/Offline/");
            ITHACAstream::read_fields(PfieldData, _p(), "./ITHACAoutput/Offline/");
            ITHACAstream::read_fields(nutFieldsData, _nut(), "./ITHACAoutput/Offline/");
            /// Compute the coefficients for train
            std::cout << "Computing the coefficients for U" << std::endl;
            Eigen::MatrixXd coeffL2U_data = ITHACAutilities::getCoeffs(UfieldData,Umodes, 0, true);
            std::cout << "Computing the coefficients for p" << std::endl;
            Eigen::MatrixXd coeffL2P_data = ITHACAutilities::getCoeffs(PfieldData,Pmodes,0, true);
            std::cout << "Computing the coefficients for nuT" << std::endl;
            Eigen::MatrixXd coeffL2Nut_data = ITHACAutilities::getCoeffs(nutFieldsData, nutModes, 0, true);
            coeffL2U_data.transposeInPlace();
            coeffL2P_data.transposeInPlace();
            coeffL2Nut_data.transposeInPlace();
            //CoeffPointDisplacement_data.transposeInPlace();
            cnpy::save(coeffL2U_data, "ITHACAoutput/NN/coeffs/coeffL2U.npy");
            cnpy::save(coeffL2P_data, "ITHACAoutput/NN/coeffs/coeffL2P.npy");
            cnpy::save(coeffL2Nut_data, "ITHACAoutput/NN/coeffs/coeffL2Nut.npy");
        }
    }
    // Function to eval the NN once the input is provided
    Eigen::MatrixXd UnsteadyNSPimpleNN::evalNet(Eigen::MatrixXd a)
    {
        //std::cout << "a before scaling \t" << a << std::endl;
        a = a.array() * scale_inp.array() + bias_inp.array() ;
        //a = (a.array() - bias_inp.array()) / scale_inp.array();
        //std::cout << "a after scaling \t" << a << std::endl;
        a.transposeInPlace();
        torch::Tensor xTensor = eigenMatrix2torchTensor(a);
        torch::Tensor out;
        std::vector<torch::jit::IValue> XTensorInp;
        XTensorInp.push_back(xTensor);
        out = netTorchscript.forward(XTensorInp).toTensor();
        //std::cout << "out \t" << out << std::endl;
        Eigen::MatrixXd g = torchTensor2eigenMatrix<double>(out);
        g.transposeInPlace();
        //std::cout << "g before rescaling \t" << g << std::endl;
        g = (g.array() - bias_out.array()) / scale_out.array();///????????/
        //g = g.array() * scale_out.array() + bias_out.array() ;
        //std::cout << "g after rescaling \t" << g << std::endl;
        return g;
    }


    Eigen::MatrixXd UnsteadyNSPimpleNN::evalLstm(Eigen::MatrixXd a)
    {

        torch::Tensor x_reshaped, predictions;
        //a = a.array() * scale_inp.array() + bias_inp.array() ;
        std::cout << "a before scaling \t" << a << std::endl;
        a = a.array() * scale_inp.array() + bias_inp.array() ;
        //a = (a.array() - bias_inp.array()) / scale_inp.array();
        std::cout << "a after scaling \t" << a << std::endl;
        //predictions.torch::Tensor::to(torch::kLong);
        a.transposeInPlace();
        torch::Tensor x = eigenMatrix2torchTensor(a);
        std::cout << "a To Tensor: " << x << std::endl;
        x_reshaped = torch::reshape(x, {1, 1, a.cols()});
        x_reshaped.to(torch::kLong);
        //torch::Tensor s, predictions;
        std::vector<torch::jit::IValue> XTensorInp;
        XTensorInp.push_back(x_reshaped);
        predictions = lstmTorchscript.forward(XTensorInp).toTensor();
        //y_reshaped = torch::reshape(predictions, {1, 1, a.cols()});
        torch::Tensor out = predictions.squeeze(0); 
        std::cout << "predictions" <<  predictions << std::endl;
        Eigen::MatrixXd g = torchTensor2eigenMatrix<double>(predictions);
        g.transposeInPlace();
        std::cout << "g before rescaling \t" << g << std::endl;
        g = (g.array() - bias_out.array()) / scale_out.array();
        std::cout << "g after rescaling: " << g << std::endl;

        g.array();

        return g;
    }

    void UnsteadyNSPimpleNN::truthSolve(List<scalar> mu_now, fileName folder = "./ITHACAoutput/Offline/")
    {

        Time& runTime = _runTime();
        dynamicFvMesh& mesh = meshPtr();
        fv::options& fvOptions = _fvOptions();
        pimpleControl& pimple = _pimple();
        volScalarField& p = _p();
        volVectorField& U = _U();
        surfaceScalarField& phi = _phi();
        IOMRFZoneList& MRF = _MRF();
        singlePhaseTransportModel& laminarTransport = _laminarTransport();
        autoPtr<incompressible::turbulenceModel> turbulence;
        turbulence = autoPtr<incompressible::turbulenceModel>
        (
             incompressible::turbulenceModel::New(U, phi, laminarTransport)
        );

        volScalarField& nut  = _nut(); 

        instantList Times = runTime.times();
        runTime.setEndTime(finalTime);
        runTime.setTime(Times[1], 1);
        const scalar ramp = 1.0;
        runTime.setDeltaT(timeStep);
        nextWrite = startTime; // timeStep initialization
        dictionary dictCoeffs(dyndict->findDict("sixDoFRigidBodyMotionCoeffs"));
        sixDoFRigidBodyMotion sDRBM(dictCoeffs, dictCoeffs, runTime);
        Foam::functionObjects::forces fomforces("fomforces", mesh, dictCoeffs);
        vector rotationAngle(Zero);
     
#include "createUfIfPresent.H"

        turbulence->validate();
        Info<< "\nStarting time loop\n" << endl;
        while (runTime.run())
        {

//#include "readDyMControls.H"
#include "CourantNo.H"
//#include "setDeltaT.H"
            runTime.setEndTime(finalTime);
            ++runTime;
            // Info<< "Value = " << runTime.deltaTValue() << nl << endl;
            Info<< "Time = " << runTime.timeName() << nl << endl;
            // --- Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())
            {
                if (pimple.firstIter() || moveMeshOuterCorrectors)
                {
                    // Do any mesh changes
                    //mesh.controlledUpdate();
                    // The following line remplace the above controlledUpdate() method
                    fomforces.calcForcesMoment();
                    sDRBMS().solve();
                    rotationAngle = quaternion(sDRBMS().motion().orientation()).eulerAngles(quaternion::XYZ);
                    mesh.movePoints(sDRBMS().curPoints());
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
                }

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
            nut = turbulence->nut();
         
            if (checkWrite(runTime))
            {
                fomforcey.append(fomforces.forceEff().y());
                fomforcex.append(fomforces.forceEff().x());

                Rx.append(sDRBMS().motion().centreOfRotation().y() ); 
                centerofmassy.append(sDRBMS().motion().centreOfMass().y() ); 
                Ry.append(rotationAngle.z() );
                Rz.append(fomforces.momentEff().z());
                //Radius.append(magSqr(R) );                   
                ITHACAstream::exportSolution(U, name(counter), folder);
                ITHACAstream::exportSolution(p, name(counter), folder);
                ITHACAstream::exportSolution(nut, name(counter), folder);
                ITHACAstream::exportSolution(phi, name(counter), folder);
                ITHACAstream::exportSolution(sDRBMS().pointDisplacement(), name(counter), folder);
                ITHACAstream::writePoints(meshPtr().points(), folder, name(counter) + "/polyMesh/");

                std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
                Ufield.append(U.clone());
                Pfield.append(p.clone());
                nutFields.append(turbulence->nut() );
                phiFields.append(phi.clone());
                Dfield.append(sDRBMS().pointDisplacement().clone());
                counter++;
                nextWrite += writeEvery;

                writeMu(mu_now);
                // --- Fill in the mu_samples with parameters (time, mu) to be used for the PODI sample points
                mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size() + 1);
                mu_samples(mu_samples.rows() - 1, 0) = atof(runTime.timeName().c_str());

                for (label i = 0; i < mu_now.size(); i++)
                {
                    mu_samples(mu_samples.rows() - 1, i + 1) = mu_now[i];
                }
            }
         

        }
        // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
        if (mu.cols() == 0)
        {
            mu.resize(1, 1);
        }

        if (mu_samples.rows() == mu.cols())
        {
            ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",folder);
        }

    } 