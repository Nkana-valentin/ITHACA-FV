// Solve the Momentum equation
// MRF.correctBoundaryVelocity(U);

// tmp<fvVectorMatrix> tUEqn
// (
//    fvm::ddt(U) 
//   -fvm::laplacian(nu, U) 
//   -fvc::div(nu * dev2(T(fvc::grad(U))))
//  ==
//    -fvc::div(phi, U)
// );
// fvVectorMatrix& UEqn = tUEqn.ref();

// UEqn.relax();

// fvOptions.constrain(UEqn);

// if (pimple.momentumPredictor())
// {
//     solve(UEqn == -fvc::grad(p) );
//     fvOptions.correct(U);
// }
// still ok
/////////////////////////////////////////////////////////////
fvVectorMatrix UEqn1
(
   fvm::ddt(U) + turbulence->divDevReff(U) 
   ==  
   T + fvOptions(U)
);
//fvVectorMatrix& UEqn = tUEqn.ref();

UEqn1.relax();
fvOptions.constrain(UEqn1);

if (pimple.momentumPredictor())
{
    solve(UEqn1 == -fvc::grad(p) );
    //UEqn1.solve();
    fvOptions.correct(U);
}


