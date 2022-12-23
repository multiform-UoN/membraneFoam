/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "membranePressureGradFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::membranePressureGradFvPatchField::membranePressureGradFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchField<scalar>(p, iF),
    K0_(p.size(), scalar(1e-9)),
    phi0_(p.size(), scalar(0.5)),
    L_(p.size(), scalar(1))
{}


Foam::membranePressureGradFvPatchField::membranePressureGradFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<scalar>(p, iF),
    K0_("membranePermeability", dict, p.size()),
    phi0_("membranePorosity", dict, p.size()),
    L_("membraneWidth", dict, p.size())
{
    this->evaluate();
}


Foam::membranePressureGradFvPatchField::membranePressureGradFvPatchField
(
    const membranePressureGradFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<scalar>(ptf, p, iF, mapper),
    K0_(mapper(ptf.K0_)),
    phi0_(mapper(ptf.phi0_)),
    L_(mapper(ptf.L_))
{}


Foam::membranePressureGradFvPatchField::membranePressureGradFvPatchField
(
    const membranePressureGradFvPatchField& ptf
)
:
    fixedGradientFvPatchField<scalar>(ptf),
    K0_(ptf.K0_),
    phi0_(ptf.phi0_),
    L_(ptf.L_)
{}


Foam::membranePressureGradFvPatchField::membranePressureGradFvPatchField
(
    const membranePressureGradFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchField<scalar>(ptf, iF),
    K0_(ptf.K0_),
    phi0_(ptf.phi0_),
    L_(ptf.L_)
{
    this->evaluate();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::membranePressureGradFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvMesh& mesh = this->internalField().mesh();
    if ( mesh.objectRegistry::template foundObject<surfaceScalarField>("phi") )
    {

        const fvsPatchField<scalar>& phi =
          patch().lookupPatchField<surfaceScalarField, scalar>("phi");

        gradient() =
            (
              -(phi / patch().magSf())
              /
              (K0_)
            );
    }
    else
    {
        gradient() = scalar(0) /  (K0_);
    }

    fixedGradientFvPatchField<scalar>::updateCoeffs();
}

void Foam::membranePressureGradFvPatchField::write(Ostream& os) const
{
    fixedGradientFvPatchField<scalar>::write(os);
    writeEntry(os, "membranePermeability", K0_);
    writeEntry(os, "membranePorosity", phi0_);
    writeEntry(os, "membraneWidth", L_);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //


namespace Foam
{
  makePatchTypeField
  (
    fvPatchScalarField,
    membranePressureGradFvPatchField
  );
}


// ************************************************************************* //
