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
    Foam::RobinPhiFvPatchScalarField

Description
    This boundary condition provide a Robin BC with an additional flux read from a surfaceVectorField

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    \endverbatim


SourceFiles
    RobinPhiFvPatchScalarField.C

Authors:
    Matteo Icardi, Federico Municchi, Nottingham (2019)
\*---------------------------------------------------------------------------*/

#ifndef RobinPhiFvPatchScalarField_H
#define RobinPhiFvPatchScalarField_H

#include "RobinFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class RobinPhiFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class RobinPhiFvPatchScalarField
:
    public RobinFvPatchScalarField
{
private:

    //- Flux variable
    word phiName_;

    //- Effective RobinK_
    scalarField RobinKeff_;

public:

    //- Runtime type information
    TypeName("RobinPhi");


    // Constructors

        //- Construct from patch and internal field
        RobinPhiFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        RobinPhiFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given RobinPhiFvPatchScalarField
        //  onto a new patch
        RobinPhiFvPatchScalarField
        (
            const RobinPhiFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        RobinPhiFvPatchScalarField
        (
            const RobinPhiFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new RobinPhiFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        RobinPhiFvPatchScalarField
        (
            const RobinPhiFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new RobinPhiFvPatchScalarField(*this, iF)
            );
        }


        // Mapping functions

             //- Map (and resize as needed) from self given a mapping object
             virtual void autoMap
             (
                 const fvPatchFieldMapper&
             );

             //- Reverse map the given fvPatchField onto this fvPatchField
             virtual void rmap
             (
                 const fvPatchField<scalar>&,
                 const labelList&
             );


        // Member functions

            virtual scalarField& RobinK()
            {
                return RobinKeff_;
            }

            virtual const scalarField& RobinK() const
            {
                return RobinKeff_;
            }

            // //- Evaluate the patch field
            // virtual void evaluate
            // (
            //     const Pstream::commsTypes commsType =
            //         Pstream::commsTypes::blocking
            // );

            virtual void updateCoeffs();

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
