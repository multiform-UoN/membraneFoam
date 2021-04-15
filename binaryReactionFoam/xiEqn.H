// Concentration Evolution Equation
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix xiEqn
    (
        fvm::ddt( ( scalar(1) - psi ) , xi )
        + fvm::div( (phi * (scalar(1) - psif) ) , xi , "div(phi,xi)")
        - fvm::div( D * fvc::snGrad(psi) * mesh.magSf() , xi , "div(phi,xi)")
        - fvm::laplacian( D*(scalar(1)-psi) , xi)
    );

    xiEqn.relax();
    xiEqn.solve();

    xi.relax();
    xi.correctBoundaryConditions();
}
