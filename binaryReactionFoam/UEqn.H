MRF.correctBoundaryVelocity(U);  //added

// Momentum predictor
tmp<fvVectorMatrix> tUEqn
(
   // Stokes approximation
   fvm::ddt(U)
   +
   fvm::div(phi, U)
   -
   fvm::laplacian(nu, U)
   -
   fvc::div(nu*dev2(fvc::grad(U)().T()))
   ==
   fvOptions(U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();
fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}
