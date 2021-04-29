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
