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

#include "membranePressureFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::membranePressureFvPatchField::membranePressureFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    K0_(p.size(), scalar(1e-9)),
    phi0_(p.size(), scalar(0.5)),
    p0_(p.size(), scalar(1))
{}


Foam::membranePressureFvPatchField::membranePressureFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    K0_("membranePermeability", dict, p.size()),
    phi0_("membranePorosity", dict, p.size()),
    p0_("outsidePressure", dict, p.size())
{
    this->evaluate();
}


Foam::membranePressureFvPatchField::membranePressureFvPatchField
(
    const membranePressureFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    K0_(mapper(ptf.K0_)),
    phi0_(mapper(ptf.phi0_)),
    p0_(mapper(ptf.p0_))
{}


Foam::membranePressureFvPatchField::membranePressureFvPatchField
(
    const membranePressureFvPatchField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    K0_(ptf.K0_),
    phi0_(ptf.phi0_),
    p0_(ptf.p0_)
{}


Foam::membranePressureFvPatchField::membranePressureFvPatchField
(
    const membranePressureFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    K0_(ptf.K0_),
    phi0_(ptf.phi0_),
    p0_(ptf.p0_)
{
    this->evaluate();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::membranePressureFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    scalarField pp( p0_ + scalar(0) / K0_ );

    const fvMesh& mesh = this->internalField().mesh();
    if ( mesh.objectRegistry::template foundObject<surfaceScalarField>("phi") )
    {

        const fvsPatchField<scalar>& phi =
          patch().lookupPatchField<surfaceScalarField, scalar>("phi");

        pp =
          (
            Foam::max
            (
              p0_
              +
              (
                  (
                    ( phi / patch().magSf() )
                  )
                /
                (K0_)
              )
            ,
            p0_
          )
        );

    }

    operator==(pp);

    fixedValueFvPatchField<scalar>::updateCoeffs();
}

void Foam::membranePressureFvPatchField::write(Ostream& os) const
{
    fixedValueFvPatchField<scalar>::write(os);
    writeEntry(os, "membranePermeability", K0_);
    writeEntry(os, "membranePorosity", phi0_);
    writeEntry(os, "outsidePressure", p0_);
}


// ************************************************************************* //


namespace Foam
{
  makePatchTypeField
  (
    fvPatchScalarField,
    membranePressureFvPatchField
  );
}


// ************************************************************************* //
