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
    Foam::roughnessModels::rModel2

Description

SourceFiles
    rModel2.C

\*---------------------------------------------------------------------------*/

#ifndef rModel2_H
#define rModel2_H

#include "roughnessModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace roughnessModels
{

/*---------------------------------------------------------------------------*\
                           Class rModel2 Declaration
\*---------------------------------------------------------------------------*/

class rModel2
:
    public roughnessModel
{
private:

    // Private data

        //- Residual Reynolds Number
        dimensionedScalar Cd_;
        dimensionedScalar k_;
        dimensionedScalar AfbyAt_;
        volScalarField alpha_;


public:

    //- Runtime type information
    TypeName("rModel2");


    // Constructors

        //- Construct from a dictionary and a phase pair
        rModel2
        (
            const volVectorField& U,
            const dictionary& roughnessModelDict
        );


    //- Destructor
    virtual ~rModel2();


    // Member Functions

        //- Drag coefficient
        virtual tmp<fvVectorMatrix> force(volVectorField& U) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace roughnessModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
