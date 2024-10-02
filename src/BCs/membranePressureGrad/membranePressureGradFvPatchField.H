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
    Foam::membranePressureGradFvPatchField

Description
    This boundary condition provides a fixed gradient condition for a membrane.

Usage
    \table
        Property     | Description             | Required    | Default value
        uniformGradient | uniform gradient     | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            membranePressureGrad;
        membranePermeability 1e-9;
        membranePorosity 1e-9;
        membraneWidth 1e-9;
    }
    \endverbatim

See also
    Foam::Function1Types
    Foam::fixedGradientFvPatchField

SourceFiles
    membranePressureGradFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef membranePressureGradFvPatchField_H
#define membranePressureGradFvPatchField_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class membranePressureGradFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class membranePressureGradFvPatchField
:
    public fixedGradientFvPatchField<scalar>
{
    // Private Data

        //- permeability
        Field<scalar> K0_;
        //- porosity
        Field<scalar> phi0_;
        //- width
        Field<scalar> L_;


public:

    //- Runtime type information
    TypeName("membranePressureGrad");


    // Constructors

        //- Construct from patch and internal field
        membranePressureGradFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        membranePressureGradFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given membranePressureGradFvPatchField
        //  onto a new patch
        membranePressureGradFvPatchField
        (
            const membranePressureGradFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        membranePressureGradFvPatchField
        (
            const membranePressureGradFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new membranePressureGradFvPatchField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        membranePressureGradFvPatchField
        (
            const membranePressureGradFvPatchField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new membranePressureGradFvPatchField(*this, iF)
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
