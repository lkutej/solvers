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
    Foam::pPFPorosityModels::myGeo

Description

SourceFiles
    myGeo.C

\*---------------------------------------------------------------------------*/

#ifndef myGeo_H
#define myGeo_H

#include "pPFPorosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pPFPorosityModels
{

/*---------------------------------------------------------------------------*\
                           Class myGeo Declaration
\*---------------------------------------------------------------------------*/

class myGeo
:
    public pPFPorosityModel
{
private:

    // Private data

        dimensionedScalar Cx_;
        dimensionedScalar Cy_;
        dimensionedScalar Cz_;
        dimensionedScalar nu_;
        dimensionedScalar C1_;
        dimensionedScalar C2_;
        dimensionedScalar DmMode_;

        volScalarField Dm_;
        volVectorField fi_;

public:

    //- Runtime type information
    TypeName("myGeo");


    // Constructors

        //- Construct from a dictionary and a phase pair
        myGeo
        (
            const volVectorField& U,
            const volScalarField& alpha,
            const dictionary& pPFPorosityModelDict
        );


    //- Destructor
    virtual ~myGeo();


    // Member Functions

        //- Drag coefficient
        virtual tmp<fvVectorMatrix> force(volVectorField& U);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pPFPorosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
