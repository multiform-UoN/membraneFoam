// Concentration Evolution Equation
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix chiEqn // check if it is correct to use volfraction as below
    (
         fvm::ddt( ( scalar(1) - psi ) , chi)
         + fvm::div( (phi * (scalar(1)-psif) ) , chi , "div(phi,chi)")
         //- fvm::div( D * fvc::snGrad(psi) * mesh.magSf() , chi , "div(phi,chi)")
         - fvm::laplacian( D * (scalar(1)-psi) , chi)
         ==
         - fvm::Sp(scalar(2)*chi*(scalar(1)-psi)*Kf,chi)
         + Kf*chi*chi*(scalar(1)-psi)
         + Kf*xi*xi*(scalar(1)-psi)
    );

    chiEqn.relax();
    chiEqn.solve();
}

/*
we have used here a linearisation of the bimolecular reaction term as follows

a*b = chi^2 - xi^2 = chi0^2 + (2chi0) (chi-chi0) - xi^2 

Rimp = -2chi0 
Rexp = -chi0^2 + 2chi0^2 + xi^2 = +chi0^2 + xi^2 = (chi0^2+xi^2) =  -(chi0^2-xi^2) + 2chi0^2 = - a*b - Rimp*chi0

xi = (a-b)/2
psi = (a+b)/2
xi^2-psi^2 = (-xi+psi)(xi+psi) = b*a
*/