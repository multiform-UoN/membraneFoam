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
    Foam::membranePressureFvPatchField

Description
    This boundary condition provides a fixed value condition for a membrane.

Usage
    \verbatim
    <patchName>
    {
        type            membranePressure;
        membranePermeability 1e-9;
        membranePorosity 1e-9;
        membraneWidth 1e-9;
        outsidePressure 0;
    }
    \endverbatim

SourceFiles
    membranePressureFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef membranePressureFvPatchField_H
#define membranePressureFvPatchField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class membranePressureFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class membranePressureFvPatchField
:
    public fixedValueFvPatchField<scalar>
{
    // Private Data

        //- permeability divided by viscosity and length
        Field<scalar> K0_;
        //- porosity
        Field<scalar> phi0_;
        //- outside Pressure
        Field<scalar> p0_;


public:

    //- Runtime type information
    TypeName("membranePressure");


    // Constructors

        //- Construct from patch and internal field
        membranePressureFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        membranePressureFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given membranePressureFvPatchField
        //  onto a new patch
        membranePressureFvPatchField
        (
            const membranePressureFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        membranePressureFvPatchField
        (
            const membranePressureFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new membranePressureFvPatchField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        membranePressureFvPatchField
        (
            const membranePressureFvPatchField&,
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
                new membranePressureFvPatchField(*this, iF)
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
