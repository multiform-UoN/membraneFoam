// Concentration Evolution Equation
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix chiEqn
    (
         fvm::ddt( ( scalar(1) - psi ) , chi)
         + fvm::div( (phi * (scalar(1)-psif) ) , chi , "div(phi,chi)")
         //- fvm::div( D * fvc::snGrad(psi) * mesh.magSf() , chi , "div(phi,chi)")
         - fvm::laplacian( D * (scalar(1)-psi) , chi)
         ==
         - fvm::Sp(scalar(2)*(scalar(1)-psi)*Kf,chi)
         + Kf*chi*chi*(scalar(1)-psi)
    );

    chiEqn.relax();
    chiEqn.solve();

    chi.relax();
    chi.correctBoundaryConditions();
}
