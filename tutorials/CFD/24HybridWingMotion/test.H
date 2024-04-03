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

romforces.calcForcesMoment();
//problem->sDRBMS().solve();
/// Solving the sixRigidMotion problem
sDRBM.update
(
     firstIter,
     ramp*(romforces.forceEff() + sDRBM.mass()*g.value()),
     ramp
    *(
        romforces.momentEff() + sDRBM.mass()*(sDRBM.momentArm() ^ g.value())
     ),
     t.deltaTValue(),
     t.deltaT0Value()
);

rotationAngle = quaternion(sDRBM.orientation()).eulerAngles(quaternion::XYZ);

std::vector<double> s(2);
s[0] = sDRBM.centreOfRotation().y(); //sDRBM.centreOfMass().y(); 
s[1] = rotationAngle.z();
for (int i = 0; i < NmodesDproj; i++)
{
    //Online value of the parameter
    pdCoeffNew(i) = problem->rbfSplines[i]->eval(s);
}

Dmodes.reconstruct(pointDisplacement, pdCoeffNew, "pointDisplacement");
problem->sDRBMS().pointDisplacement().primitiveFieldRef() = pointDisplacement.primitiveFieldRef();
pointConstraints::New
(
    problem->sDRBMS().pointDisplacement().mesh()

).constrainDisplacement(problem->sDRBMS().pointDisplacement() );
//// Move the new current points
mesh.movePoints(problem->sDRBMS().curPoints());                            