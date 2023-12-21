void solve()
 {
     const Time& t = mesh().time();
 
     if (mesh().nPoints() != points0().size())
     {
         FatalErrorInFunction
             << "The number of points in the mesh seems to have changed." << endl
             << "In constant/polyMesh there are " << points0().size()
             << " points; in the current mesh there are " << mesh().nPoints()
             << " points." << exit(FatalError);
     }
 
     // Store the motion state at the beginning of the time-stepbool
     bool firstIter = false;
     if (curTimeIndex_ != t.timeIndex())
     {
         newTime();
         curTimeIndex_ = t.timeIndex();
         firstIter = true;
     }
 
     dimensionedVector g("g", dimAcceleration, Zero);
 
     if (mesh().foundObject<uniformDimensionedVectorField>("g"))
     {
         g = mesh().lookupObject<uniformDimensionedVectorField>("g");
     }
     else if (coeffDict().found("g"))
     {
         coeffDict().lookup("g") >> g;
     }
 
     // scalar ramp = min(max((t.value() - 5)/10, 0), 1);
     scalar ramp = 1.0;
 
     if (test_)
     {
         update
         (
             firstIter,
             ramp*(mass()*g.value()),
             ramp*(mass()*(momentArm() ^ g.value())),
             t.deltaTValue(),
             t.deltaT0Value()
         );
     }
     else
     {
         dictionary forcesDict;
 
         forcesDict.add("type", functionObjects::forces::typeName);
         forcesDict.add("patches", patches_);
         forcesDict.add("rhoInf", rhoInf_);
         forcesDict.add("rho", rhoName_);
         forcesDict.add("CofR", centreOfRotation());
 
         functionObjects::forces f("forces", t, forcesDict);
 
         f.calcForcesMoment();
 
         update
         (
             firstIter,
             ramp*(f.forceEff() + mass()*g.value()),
             ramp
            *(
                f.momentEff()
              + mass()*(momentArm() ^ g.value())
             ),
             t.deltaTValue(),
             t.deltaT0Value()
         );
     }
 
     // Update the displacements
     pointDisplacement_.primitiveFieldRef() =
         transform(points0(), scale_) - points0();
 
     // Displacement has changed. Update boundary conditions
     pointConstraints::New
     (
         pointDisplacement_.mesh()
     ).constrainDisplacement(pointDisplacement_);
 }