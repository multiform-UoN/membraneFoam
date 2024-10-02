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
    Foam::membraneVelocityFvPatchField

Description
    This boundary condition provides a fixed value condition for a membrane.

Usage
    \verbatim
    <patchName>
    {
        type            membraneVelocity;
        membranePermeability 1e-9;
        membranePorosity 1e-9;
        membraneWidth 1e-9;
        outsidePressure 0;
    }
    \endverbatim

SourceFiles
    membraneVelocityFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef membraneVelocityFvPatchField_H
#define membraneVelocityFvPatchField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class membraneVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class membraneVelocityFvPatchField
:
    public fixedValueFvPatchField<vector>
{
    // Private Data

        //- initial permeability
        scalar K0_;
        //- initial porosity
        scalar phi0_;
        //- outside (permeate) pressure
        Field<scalar> p0_;
        //- width
        scalar w_;
        //- permeability
        scalarField K_;
        //- scaling permeability
        scalar Ks_;
        //- porosity
        scalarField phi_;

        //- permeability models
        word solidM_;   // model
        word solidS_;   // solid name
        scalar solidV_; // molar volume
        bool binaryReaction_; // binaryReaction BC

        //- osmotic pressure
        word osmoticP_;   // osmotic pressure species name
        scalar osmoticC_; // osmotic coefficient

        // - debug and plotting
        bool writeAvg_;




public:

    //- Runtime type information
    TypeName("membraneVelocity");


    // Constructors

        //- Construct from patch and internal field
        membraneVelocityFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        membraneVelocityFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given membraneVelocityFvPatchField
        //  onto a new patch
        membraneVelocityFvPatchField
        (
            const membraneVelocityFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        membraneVelocityFvPatchField
        (
            const membraneVelocityFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return tmp<fvPatchField<vector>>
            (
                new membraneVelocityFvPatchField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        membraneVelocityFvPatchField
        (
            const membraneVelocityFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector>>
            (
                new membraneVelocityFvPatchField(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#endif

// ************************************************************************* //
