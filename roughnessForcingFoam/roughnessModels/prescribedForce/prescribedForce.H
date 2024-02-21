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
    Foam::roughnessModels::prescribedForce

Description

SourceFiles
    prescribedForce.C

\*---------------------------------------------------------------------------*/

#ifndef prescribedForce_H
#define prescribedForce_H

#include "roughnessModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace roughnessModels
{

/*---------------------------------------------------------------------------*\
                           Class prescribedForce Declaration
\*---------------------------------------------------------------------------*/

class prescribedForce
:
    public roughnessModel
{
private:

    // Private data

        //- Residual Reynolds Number
        volScalarField dragCoeff_;

public:

    //- Runtime type information
    TypeName("prescribedForce");


    // Constructors

        //- Construct from a dictionary and a phase pair
        prescribedForce
        (
            const volVectorField& U,
            const dictionary& roughnessModelDict
        );


    //- Destructor
    virtual ~prescribedForce();


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
