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
