// Concentration Evolution Equation
while (pimple.correctNonOrthogonal())
{
      fvScalarMatrix psiEqn
    (

         fvm::ddt(psi)
         + fvm::div(phic,psi)
         - fvm::laplacian(Dpsi,psi)
         ==
         - fvm::Sp(Kf*chi*chi/psifact,psi)
         + Kf*chi*chi/psifact
    );

    psiEqn.relax();
    psiEqn.solve();

    psi.relax();
    psi.correctBoundaryConditions();
}
