// Concentration Evolution Equation
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix chiEqn
    (
         fvm::ddt( chi)
         + fvm::div( (phi ) , chi , "div(phi,chi)")
         - fvm::laplacian( D  , chi)
    );

    chiEqn.relax();
    chiEqn.solve();

    fvScalarMatrix chiEqn2
    (
         fvm::ddt( chi2)
         + fvm::div( (phi ) , chi2 , "div(phi,chi)")
         //- fvm::div( D * fvc::snGrad(psi) * mesh.magSf() , chi , "div(phi,chi)")
         - fvm::laplacian( D , chi2)
    );

    chiEqn2.relax();
    chiEqn2.solve();

}
