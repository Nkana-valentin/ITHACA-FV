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

tmp<fvVectorMatrix> tUEqn
(
   fvm::ddt(U) 
  //-fvm::laplacian(nu, U) 
  //-fvc::div(nu * dev2(T(fvc::grad(U))))
 ==
   -fvc::div(phi, U) + fvc::div(nu, mesh.Cf()) + nu * dev2(T(fvc::grad(U)))
);
fvVectorMatrix& UEqn = tUEqn.ref();

//UEqn.relax();

//fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p) );
    //UEqn.solve();
    //fvOptions.correct(U);
}


