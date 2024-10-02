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



Class
    Foam::noFluxFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a base class for 'noFlux' type boundary
    conditions setting all coefficients to zero

SourceFiles
    noFluxFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef noFluxFvPatchField_H
#define noFluxFvPatchField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class noFluxFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class noFluxFvPatchField
:
    public fvPatchField<Type>
{
    // Private data

        //- Source term
        //Field<Type> noFluxF_;

public:

    //- Runtime type information
    TypeName("noFlux");


    // Constructors

        //- Construct from patch and internal field
        noFluxFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        noFluxFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given noFluxFvPatchField onto a new patch
        noFluxFvPatchField
        (
            const noFluxFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );

        //- Construct as copy
        noFluxFvPatchField
        (
            const noFluxFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new noFluxFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        noFluxFvPatchField
        (
            const noFluxFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new noFluxFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Return false: this patch field is not altered by assignment
            virtual bool assignable() const
            {
                return false;
            }


        // Return defining fields

            // virtual Field<Type>& noFluxF()
            // {
            //     return noFluxF_;
            // }
            //
            // virtual const Field<Type>& noFluxF() const
            // {
            //     return noFluxF_;
            // }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const;

            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const UList<Type>&) {}

        virtual void operator=(const fvPatchField<Type>&) {}
        virtual void operator+=(const fvPatchField<Type>&) {}
        virtual void operator-=(const fvPatchField<Type>&) {}
        virtual void operator*=(const fvPatchField<scalar>&) {}
        virtual void operator/=(const fvPatchField<scalar>&) {}

        virtual void operator+=(const Field<Type>&) {}
        virtual void operator-=(const Field<Type>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const Type&) {}
        virtual void operator+=(const Type&) {}
        virtual void operator-=(const Type&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "noFluxFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
