// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) 
    + fvm::div(phi, U)
    + MRF.DDt(U)
    + turbulence->divDevReff(U)
    //+ fvc::grad(p) 
    == -fvc::grad(p) + fvOptions(U) 
);
fvVectorMatrix& UEqn = tUEqn.ref();
vectorField resF(UEqn.residual());

for (label i = 0; i < R.internalField().size(); i++)
{
    R.ref()[i] = resF[i];
}

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    //solve(UEqn == -fvc::grad(p));
    UEqn.solve();

    fvOptions.correct(U);
}
