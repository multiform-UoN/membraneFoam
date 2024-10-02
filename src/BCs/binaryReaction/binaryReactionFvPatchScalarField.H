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
    Foam::binaryReactionFvPatchScalarField

Description
    This boundary condition provide a reversible binary reactive flux for a the
    reactive component with solid formed on the surface

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    \endverbatim


SourceFiles
    binaryReactionFvPatchScalarField.C

Authors:
    Matteo Icardi, Nottingham (2020)
\*---------------------------------------------------------------------------*/

#ifndef binaryReactionFvPatchScalarField_H
#define binaryReactionFvPatchScalarField_H

#include "RobinPhiFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class binaryReactionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class binaryReactionFvPatchScalarField
:
    public RobinPhiFvPatchScalarField
{
private:

    //- Boundary concentration
    scalarField S_;

    //- Boundary concentration
    //  at previous time step
    scalarField S0_;

    //- Effective RobinK_
    scalarField RobinKeff_;

    //- Reversible reaction rate Kd
    scalarField Kd_;

    //- Effective RobinF_
    scalarField RobinFeff_;

    //- Current time index.
    //  Used to trigger the storing of the old-time value
    mutable label timeIndex_;

    //- non reactive component
    word xiName_;

    word reactionType_;

    // - debug and plotting
    bool writeAvg_;

public:

    //- Runtime type information
    TypeName("binaryReaction");


    // Constructors

        //- Construct from patch and internal field
        binaryReactionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        binaryReactionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given binaryReactionFvPatchScalarField
        //  onto a new patch
        binaryReactionFvPatchScalarField
        (
            const binaryReactionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        binaryReactionFvPatchScalarField
        (
            const binaryReactionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new binaryReactionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        binaryReactionFvPatchScalarField
        (
            const binaryReactionFvPatchScalarField&,
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
                new binaryReactionFvPatchScalarField(*this, iF)
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

            scalar Savg() const
            {
                return gSum(this->patch().magSf()*S_)/gSum(this->patch().magSf());
            }

            scalarField S() const
            {
                return S_;
            }

            void rescaleS(scalar avg)
            {
                S_ = S_*avg/Savg();
            }

            virtual scalarField& RobinF()
            {
                return RobinFeff_;
            }

            virtual const scalarField& RobinF() const
            {
                return RobinFeff_;
            }

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

            //- Evaluate the patch field
            virtual void updateCoeffs();

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
