/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
\\    /   O peration     |
\\  /    A nd           | Copyright (C) 2015-2019
\\/     M anipulation  | Matteo Icardi, Federico Municchi
-------------------------------------------------------------------------------
License
This file is derivative work of OpenFOAM.

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

#include "binaryReactionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

#define STEADY_TIMESTEP 1e10

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binaryReactionFvPatchScalarField::binaryReactionFvPatchScalarField
(
  const fvPatch& p,
  const DimensionedField<scalar, volMesh>& iF
)
:
RobinPhiFvPatchScalarField(p, iF),
S_(p.size()),
S0_(p.size()),
RobinKeff_(p.size()),
Kd_(p.size()),
RobinFeff_(p.size()),
timeIndex_(-1),
xiName_("xi"),
reactionType_("binary")
{

}


Foam::binaryReactionFvPatchScalarField::binaryReactionFvPatchScalarField
(
  const fvPatch& p,
  const DimensionedField<scalar, volMesh>& iF,
  const dictionary& dict
)
:
RobinPhiFvPatchScalarField(p, iF, dict),
S_("S",dict,p.size()),
S0_(S_),
RobinKeff_(p.size(),scalar(0)),
Kd_("Kd",dict,p.size()),
RobinFeff_(p.size(),scalar(0)),
timeIndex_(-1),
xiName_(dict.lookupOrDefault<word>("xi", "xi")),
reactionType_(dict.lookupOrDefault<word>("reactionType", "binary"))
{
}


Foam::binaryReactionFvPatchScalarField::binaryReactionFvPatchScalarField
(
  const binaryReactionFvPatchScalarField& ptf,
  const fvPatch& p,
  const DimensionedField<scalar, volMesh>& iF,
  const fvPatchFieldMapper& mapper
)
:
RobinPhiFvPatchScalarField(p, iF),
S_(mapper(ptf.S_)),//,mapper),
S0_(mapper(ptf.S0_)),//,mapper),
RobinKeff_(mapper(ptf.RobinKeff_)),//,mapper),
Kd_(mapper(ptf.Kd_)),//,mapper),
RobinFeff_(mapper(ptf.RobinFeff_)),//,mapper),
timeIndex_(ptf.timeIndex_),
xiName_(ptf.xiName_),
reactionType_(ptf.reactionType_)
{

}


Foam::binaryReactionFvPatchScalarField::binaryReactionFvPatchScalarField
(
  const binaryReactionFvPatchScalarField& ptf
)
:
RobinPhiFvPatchScalarField(ptf),
S_(ptf.S_),
S0_(ptf.S0_),
RobinKeff_(ptf.RobinKeff_),
Kd_(ptf.Kd_),
RobinFeff_(ptf.RobinFeff_),
timeIndex_(ptf.timeIndex_),
xiName_(ptf.xiName_),
reactionType_(ptf.reactionType_)
{

}


Foam::binaryReactionFvPatchScalarField::binaryReactionFvPatchScalarField
(
  const binaryReactionFvPatchScalarField& ptf,
  const DimensionedField<scalar, volMesh>& iF
)
:
RobinPhiFvPatchScalarField(ptf, iF),
S_(ptf.S_),
S0_(ptf.S0_),
RobinKeff_(ptf.RobinKeff_),
Kd_(ptf.Kd_),
RobinFeff_(ptf.RobinFeff_),
timeIndex_(ptf.timeIndex_),
xiName_(ptf.xiName_),
reactionType_(ptf.reactionType_)
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::binaryReactionFvPatchScalarField::autoMap
(
  const fvPatchFieldMapper&   m
)
{
  RobinPhiFvPatchScalarField::autoMap(m);
  m(S_,S_);//S_.autoMap(mapper);
  m(S0_,S0_);//S0_.autoMap(mapper);
  m(RobinKeff_,RobinKeff_);//.autoMap(mapper);
  m(Kd_,Kd_);//.autoMap(mapper);
  m(RobinFeff_,RobinFeff_);//.autoMap(mapper);
}

void Foam::binaryReactionFvPatchScalarField::rmap
(
  const fvPatchField<scalar>& ptf,
  const labelList& addr
)
{
  RobinPhiFvPatchScalarField::rmap(ptf,addr);
  const binaryReactionFvPatchScalarField& mptf =
  refCast<const binaryReactionFvPatchScalarField>(ptf);
  S_.rmap(mptf.S_,addr);
  S0_.rmap(mptf.S0_,addr);
  RobinKeff_.rmap(mptf.RobinKeff_,addr);
  Kd_.rmap(mptf.Kd_,addr);
  RobinFeff_.rmap(mptf.RobinFeff_,addr);
}

void Foam::binaryReactionFvPatchScalarField::write(Ostream& os) const
{
  RobinPhiFvPatchScalarField::write(os);
  writeEntry(os, "Kd", Kd_);
  writeEntry(os, "S", S_);
  writeEntry(os, "RobinKeff2", RobinKeff_);
  writeEntry(os, "RobinFeff", RobinFeff_);
  writeEntry(os, "xi", xiName_);
  writeEntry(os, "reactionType", reactionType_);
}

void Foam::binaryReactionFvPatchScalarField::updateCoeffs()
// void Foam::binaryReactionFvPatchScalarField::evaluate
// (
//   const Pstream::commsTypes commsType
// )
{

  RobinPhiFvPatchScalarField::updateCoeffs();

  const fvMesh& mesh = this->internalField().mesh();

  word ddtScheme
  (
      mesh.ddtScheme("binaryReaction")
  );

  scalar deltaT = this->db().time().deltaTValue();

  if (ddtScheme=="steadyState")
  {
    Info << "binaryReactionBC in steady state, running with time step " << STEADY_TIMESTEP;
    deltaT = STEADY_TIMESTEP;
  }

  scalarField& C(*this);
  const scalarField& RobinKorig = RobinFvPatchScalarField::RobinK();
  const scalarField& RobinK0 = RobinPhiFvPatchScalarField::RobinK();
  const scalarField& RobinF0 = RobinPhiFvPatchScalarField::RobinF();

  if ( mesh.objectRegistry::template foundObject<volScalarField>(xiName_) )
  {
    const fvPatchField<scalar>& xi =
        patch().lookupPatchField<volScalarField, scalar>(xiName_);
  
    //- Get old time concentration
    if(timeIndex_!=this->db().time().timeIndex())
    {
      S0_ = S_;
      timeIndex_ = this->db().time().timeIndex();
    }

    scalarField reaction((C-xi)*(C+xi));  // chi^2 - xi^2

    // if (reactionType_=="2")
    // {
    //   reaction *= (C-xi)/scalar(2);
    // }

    // - Evolution equation for S (implicit time stepping)
    S_ = (S0_ + deltaT*RobinKorig*reaction )  // mol/area => RobinK0 [vol/mol l/t] = K0_bulk*width
          /
          (scalar(1) + deltaT*Kd_);

    //- Firt guess for Robin coefficients
    RobinKeff_ = - RobinKorig + RobinK0   // advective term
                 - scalar(2)*RobinKorig*C; // linearised reactive term (implicit)

    RobinFeff_ =   RobinKorig*xi*xi // K chi0^2 + Kxi^2
                 + RobinKorig*C*C
                 + Kd_*S_;
  }

  RobinFvPatchScalarField::updateCoeffs();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  makePatchTypeField
  (
    fvPatchScalarField,
    binaryReactionFvPatchScalarField
  );
}

// ************************************************************************* //
