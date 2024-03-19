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