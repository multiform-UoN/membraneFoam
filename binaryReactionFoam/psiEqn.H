// Concentration Evolution Equation
while (pimple.correctNonOrthogonal())
{
      fvScalarMatrix psiEqn
    (

         fvm::ddt(psi)
         + fvm::div(phic,psi)
         - fvm::laplacian(Dpsi,psi)
         ==
         // - fvm::Sp(Kr/psifact,psi) // TODO add reverse reaction
         Kf*chi*chi*psifact
         - Kf*xi*xi*psifact
    );

    psiEqn.relax();
    psiEqn.solve();
}

