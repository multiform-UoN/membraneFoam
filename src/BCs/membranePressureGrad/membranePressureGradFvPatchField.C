/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) The OpenFOAM Foundation, Ltd.
     \\/     M anipulation  | 
-------------------------------------------------------------------------------
  
  Modifications and additional contributions:
  
  MultiForm Group, University of Nottingham
  Copyright (C) Matteo Icardi and collaborators
  Contributions to this file are licensed under the same GPLv3 terms as OpenFOAM.

  This work is based on OpenFOAM, with substantial portions copied, modified, or 
  extended under the GNU General Public License (GPLv3).
  Please refer to the original OpenFOAM copyright notice for the base framework, 
  and to MultiForm Group for new additions or modifications. For the full terms of 
  this license, see <http://www.gnu.org/licenses/>.

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
