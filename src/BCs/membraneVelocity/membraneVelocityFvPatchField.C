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

#include "membraneVelocityFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "binaryReactionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::membraneVelocityFvPatchField::membraneVelocityFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    K0_(scalar(1e-9)),
    phi0_(scalar(0.5)),
    p0_(p.size(), scalar(0)),
    w_(scalar(1e-4)),
    K_(p.size(),K0_),
    Ks_(K0_*scalar(1e-3)),
    phi_(p.size(),scalar(0.5)),
    solidM_("none"),
    solidS_("none"),
    solidV_(scalar(1)),
    binaryReaction_(false),
    osmoticP_("none"),
    osmoticC_(scalar(1)),
    writeAvg_(false)
{}


Foam::membraneVelocityFvPatchField::membraneVelocityFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    K0_(dict.lookupOrDefault<scalar>("membranePermeability",scalar(1e-9))),
    phi0_(dict.lookupOrDefault<scalar>("membranePorosity",scalar(.5))),
    p0_("outsidePressure", dict, p.size()),
    w_(dict.lookupOrDefault<scalar>("membraneWidth",scalar(1e-4))),
    K_(p.size(), K0_),
    Ks_(dict.lookupOrDefault<scalar>("solidPermeability",scalar(1e-12))),
    phi_(p.size(), phi0_),
    solidM_(dict.lookupOrDefault<word>("solidModel","none")),
    solidS_(dict.lookupOrDefault<word>("solidSource","none")),
    solidV_(dict.lookupOrDefault<scalar>("solidVolume",scalar(1))),
    binaryReaction_(dict.lookupOrDefault<bool>("solidBinaryReaction", false)),
    osmoticP_(dict.lookupOrDefault<word>("osmoticPressure","none")),
    osmoticC_(dict.lookupOrDefault<scalar>("osmoticCoefficient",scalar(1))),
    writeAvg_(dict.lookupOrDefault<bool>("writeAvg", false))
{
    this->evaluate();
}


Foam::membraneVelocityFvPatchField::membraneVelocityFvPatchField
(
    const membraneVelocityFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    K0_((ptf.K0_)),
    phi0_((ptf.phi0_)),
    p0_(mapper(ptf.p0_)),
    w_(ptf.w_),
    K_(mapper(ptf.K_)),
    Ks_((ptf.Ks_)),
    phi_(mapper(ptf.phi_)),
    solidM_(ptf.solidM_),
    solidS_(ptf.solidS_),
    solidV_(ptf.solidV_),
    binaryReaction_(ptf.binaryReaction_),
    osmoticP_(ptf.osmoticP_),
    osmoticC_(ptf.osmoticC_),
    writeAvg_(ptf.writeAvg_)
{}


Foam::membraneVelocityFvPatchField::membraneVelocityFvPatchField
(
    const membraneVelocityFvPatchField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    K0_(ptf.K0_),
    phi0_(ptf.phi0_),
    p0_(ptf.p0_),
    w_(ptf.w_),
    K_((ptf.K_)),
    Ks_((ptf.Ks_)),
    phi_((ptf.phi_)),
    solidM_(ptf.solidM_),
    solidS_(ptf.solidS_),
    solidV_(ptf.solidV_),
    binaryReaction_(ptf.binaryReaction_),
    osmoticP_(ptf.osmoticP_),
    osmoticC_(ptf.osmoticC_),
    writeAvg_(ptf.writeAvg_)
{}


Foam::membraneVelocityFvPatchField::membraneVelocityFvPatchField
(
    const membraneVelocityFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    K0_(ptf.K0_),
    phi0_(ptf.phi0_),
    p0_(ptf.p0_),
    w_(ptf.w_),
    K_((ptf.K_)),
    Ks_((ptf.Ks_)),
    phi_((ptf.phi_)),
    solidM_(ptf.solidM_),
    solidS_(ptf.solidS_),
    solidV_(ptf.solidV_),
    binaryReaction_(ptf.binaryReaction_),
    osmoticP_(ptf.osmoticP_),
    osmoticC_(ptf.osmoticC_),
    writeAvg_(ptf.writeAvg_)
{
    this->evaluate();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::membraneVelocityFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvMesh& mesh = this->internalField().mesh();
    const scalar area(gSum(this->patch().magSf()));

    // Find the pressure
    scalarField deltaP((-p0_));
    if  (
        mesh.objectRegistry:: template foundObject<volScalarField>("p")
        )
    {
        const fvPatchField<scalar>& p =
            patch().lookupPatchField<volScalarField, scalar>("p");
        deltaP += p;
    }

    // Find the viscosity
    scalar nu(1);
    if (
        mesh.objectRegistry:: template foundObject<IOdictionary>("transportProperties")
        )
    {
        nu = dimensionedScalar(db().lookupObject<IOdictionary>
                (
                    "transportProperties"
                ).lookup("nu")).value();
    }
    else
    {
        WarningIn("membraneVelocityFvPatchField")
            << "No transportProperties found, using nu=1" << endl;
    }


    if (writeAvg_)
    {
        Info<<"membraneVelocity BC " << this->patch().name()
          << " Time = " << mesh.time().timeName()
          << " Vel = " <<  gSum(this->patch().magSf()*(*this))/area
          << endl;
    }
    if  (
        mesh.objectRegistry::
        template foundObject<volScalarField>(osmoticP_)
        )
    {
        // Concentration
        const fvPatchScalarField& cc
              =
              patch().lookupPatchField<volScalarField, scalar>(osmoticP_);
        deltaP -= osmoticC_*cc;
        // assuming zero concentration on the other side
        if (writeAvg_)
        {
            Info << " DeltaPTot = " << gSum(this->patch().magSf()*deltaP)/area        // average total pressure
                 << " DeltaPOsm = " << gSum(this->patch().magSf()*osmoticC_*cc)/area  // average osmotic pressure
                 << endl;
        }
    }

    K_ = K0_;

    if  (
        mesh.objectRegistry::
        template foundObject<volScalarField>(solidS_)
        )
    {
        //- solid concentration on the surface
        const fvPatchScalarField& cc
              =
              patch().lookupPatchField<volScalarField, scalar>(solidS_);

        scalarField ss(cc);

        // - binaryReaction BC
        if  (binaryReaction_)
        {
          const binaryReactionFvPatchScalarField& Cpr
          (
           refCast<const binaryReactionFvPatchScalarField>(cc)
          );

          ss = Cpr.S();
        }

        if (solidM_=="KozenyCarman32")
        {
            // porosity
            phi_ = phi0_ - solidV_*ss/w_;
            phi_ *= Foam::pos(phi_);

            // // Additional solid (hard-coded)
            // if  (
            //     mesh.objectRegistry::
            //     template foundObject<volScalarField>(solidS_+"2")
            //     )
            // {
            //     //- solid concentration on the surface
            //     const fvPatchScalarField& cc2
            //         =
            //         patch().lookupPatchField<volScalarField, scalar>(solidS_+"2");

            //     scalarField ss2(cc2);

            //     // - binaryReaction BC
            //     if  (binaryReaction_)
            //     {
            //         const binaryReactionFvPatchScalarField& Cpr
            //         (
            //         refCast<const binaryReactionFvPatchScalarField>(cc)
            //         );

            //         ss2 = Cpr.S();
            //     }

            //     // porosity
            //     phi_ = phi_ - solidV_*ss2/w_;
            //     phi_ *= Foam::pos(phi_);
            // }
            
            //- permeability
            K_ =  K0_
                  *
                  Foam::pow(scalar(1)-phi0_,2) * Foam::pow(phi_,3)
                  /
                  (Foam::pow(scalar(1)-phi_,2) * Foam::pow(phi0_,3));


            if (writeAvg_)
            {
                Info << " Perm = " << gSum(this->patch().magSf()*K_)/area
                    << " Por = " << gSum(this->patch().magSf()*phi_)/area << endl;
            }
        }
        else if (solidM_=="scalingLayer")
        {
            // // DEBUG
            // Info << endl << "Scaling layer model" << endl;
            // Info << "K0_ = " << K0_ << endl;
            // Info << "Ks_ = " << Ks_ << endl;
            // Info << (w_/K0_ + solidV_*ss/Ks_) << endl;


            // harmonic average of K0_ and Ks_ weighted by w_ and ss respectively
            K_ = (w_+solidV_*ss)/(w_/K0_ + solidV_*ss/Ks_); 
            // ss here takes the role of the scaling layer thickness

            if (writeAvg_)
            {
                Info << " Perm = " << gSum(this->patch().magSf()*K_)/area
                    << " scalingLayer = " << gSum(this->patch().magSf()*ss)/area << endl;
            }
        }
    }

    operator==
      (
          patch().nf()
          *
          deltaP
          *
          K_/(nu*w_)
      );

    fixedValueFvPatchField<vector>::updateCoeffs();

    if (writeAvg_)
    {
        Info << endl;
    }


}

void Foam::membraneVelocityFvPatchField::write(Ostream& os) const
{
    fixedValueFvPatchField<vector>::write(os);
    writeEntry(os, "membranePermeability", K0_);
    writeEntry(os, "membranePorosity", phi0_);
    writeEntry(os, "outsidePressure", p0_);
    writeEntry(os, "membraneWidth", w_);
    writeEntry(os, "membranePermeabilityFinal", K_);
    writeEntry(os, "membranePorosityFinal", phi_);
    writeEntry(os, "solidModel", solidM_);
    writeEntry(os, "solidSource", solidS_);
    writeEntry(os, "solidVolume", solidV_);
    writeEntry(os, "solidBinaryReaction", binaryReaction_);
    writeEntry(os, "osmoticPressure", osmoticP_);
    writeEntry(os, "osmoticCoefficient", osmoticC_);
    writeEntry(os, "writeAvg", writeAvg_);
}


// ************************************************************************* //


namespace Foam
{
  makePatchTypeField
  (
    fvPatchVectorField,
    membraneVelocityFvPatchField
  );
}


// ************************************************************************* //
