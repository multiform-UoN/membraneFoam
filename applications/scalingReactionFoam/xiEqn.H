// Concentration Evolution Equation
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix xiEqn
    (
        fvm::ddt( xi )
        + fvm::div( (phi  ) , xi , "div(phi,xi)")
        - fvm::laplacian( D, xi)
    );

    xiEqn.relax();
    xiEqn.solve();


    fvScalarMatrix xiEqn2
    (
        fvm::ddt( xi2 )
        + fvm::div( (phi  ) , xi2 , "div(phi,xi)")
        - fvm::laplacian( D, xi2)
    );

    xiEqn2.relax();
    xiEqn2.solve();
}
