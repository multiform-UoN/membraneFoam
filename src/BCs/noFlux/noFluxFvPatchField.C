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

#include "noFluxFvPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::noFluxFvPatchField<Type>::noFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF)//,
//    noFluxF_(p.size())
{
}


template<class Type>
Foam::noFluxFvPatchField<Type>::noFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict, false)//,
//    noFluxF_("noFluxF", dict, p.size())
{
    evaluate();
}


template<class Type>
Foam::noFluxFvPatchField<Type>::noFluxFvPatchField
(
    const noFluxFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    fvPatchField<Type>(ptf, p, iF, mapper, mappingRequired)//,
//    noFluxF_(mapper(ptf.noFluxF_))
{
    if (mappingRequired && notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}


template<class Type>
Foam::noFluxFvPatchField<Type>::noFluxFvPatchField
(
    const noFluxFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf)//,
//    noFluxF_(ptf.noFluxF_)
{}


template<class Type>
Foam::noFluxFvPatchField<Type>::noFluxFvPatchField
(
    const noFluxFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF)//,
//    noFluxF_(ptf.noFluxF_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::noFluxFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
//    m(noFluxF_,noFluxF_);//.autoMap(m);
}


template<class Type>
void Foam::noFluxFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const noFluxFvPatchField<Type>& mptf =
        refCast<const noFluxFvPatchField<Type>>(ptf);

//    noFluxF_.rmap(mptf.noFluxF_, addr);
}

template<class Type>
void Foam::noFluxFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        (
          this->patchInternalField()
        )
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::noFluxFvPatchField<Type>::snGrad() const
{
  return tmp<Field<Type>>
  (
      new Field<Type>(this->size(), Zero)
  );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::noFluxFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
  return tmp<Field<Type>>
  (
      new Field<Type>(this->size(), Zero)
  );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::noFluxFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
  return tmp<Field<Type>>
  (
      new Field<Type>(this->size(), Zero)
  );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::noFluxFvPatchField<Type>::gradientInternalCoeffs() const
{
  return tmp<Field<Type>>
  (
      new Field<Type>(this->size(), Zero)
  );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::noFluxFvPatchField<Type>::gradientBoundaryCoeffs() const
{
  return tmp<Field<Type>>
  (
      new Field<Type>(this->size(), Zero)
  );
}


template<class Type>
void Foam::noFluxFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
//    writeEntry(os,"noFluxF",noFluxF_);//noFluxF_.writeEntry("noFluxF", os);
    writeEntry(os,"value",*this);//this->writeEntry("value", os);
}



template<class Type>
void Foam::noFluxFvPatchField<Type>::updateCoeffs()
 {
     if (this->updated())
     {
         return;
     }

     fvPatchField<Type>::updateCoeffs();
 }

// ************************************************************************* //
