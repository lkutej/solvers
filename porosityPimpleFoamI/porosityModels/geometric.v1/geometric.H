/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::pPFPorosityModels::geometric

Description

SourceFiles
    geometric.C

\*---------------------------------------------------------------------------*/

#ifndef geometric_H
#define geometric_H

#include "pPFPorosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pPFPorosityModels
{

/*---------------------------------------------------------------------------*\
                           Class geometric Declaration
\*---------------------------------------------------------------------------*/

class geometric
:
    public pPFPorosityModel
{
private:

    // Private data

        dimensionedScalar nu_;
        dimensionedScalar C1_;
        dimensionedScalar C2_;
        dimensionedScalar N_;

        volScalarField Dm_;

public:

    //- Runtime type information
    TypeName("geometric");


    // Constructors

        //- Construct from a dictionary and a phase pair
        geometric
        (
            const volVectorField& U,
            const volScalarField& alpha,
            const dictionary& pPFPorosityModelDict
        );


    //- Destructor
    virtual ~geometric();


    // Member Functions

        //- Drag coefficient
        virtual tmp<fvVectorMatrix> force(volVectorField& U) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pPFPorosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
