volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
  + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
);

MRF.makeRelative(phiHbyA);

adjustPhi(phiHbyA, U, p);

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU(), MRF);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);
    
    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
    

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}
//List<Eigen::MatrixXd> RedLinSysP;
//bold = b;

// Non-orthogonal pressure corrector loop
// while (pimple.correctNonOrthogonal())
// {
//     fvScalarMatrix pEqn
//     (
//         fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
//     );

//     pEqn.setReference(pRefCell, pRefValue);
//     RedLinSysP = Pmodes.project(pEqn, Nmodes);
//     b = reducedProblem::solveLinearSys(RedLinSysP, b, presidual);
//     Pmodes.reconstruct(p, b, "p");

//     if (pimple.finalNonOrthogonalIter())
//     {
//         phi = phiHbyA - pEqn.flux();
//     }
// }    

//#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();
U = HbyA - rAtU()*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
// Eval the non linear term
//T = -fvc::div(phi, U);
// Online  coeffs
// thetaon = problem->DEIMObject->onlineCoeffs(T);
// // std::cout << problem->DEIMObject->MatrixOnline.cols() << std::endl;
// // std::cout << problem->DEIMObject->MatrixOnline.rows() << std::endl;
// //  std::cout << "=============================================================" << std::endl;
// // std::cout << thetaon.cols() << std::endl;
// // std::cout << thetaon.rows() << std::endl;
//  std::cout << "=============================================================" << std::endl;
// //Eigen::VectorXd c = problem->DEIMObject->MatrixOnline * thetaon;
//  Eigen::VectorXd c = problem->ReducedVectors * thetaon;
//T = Foam2Eigen::Eigen2field(T, c);
//problem->Tmodes.reconstruct(T, c, "T");


