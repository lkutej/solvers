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
    Foam::roughnessModels::prescribedDragAve

Description

SourceFiles
    prescribedDragAve.C

\*---------------------------------------------------------------------------*/

#ifndef prescribedDragAve_H
#define prescribedDragAve_H

#include "roughnessModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace roughnessModels
{

/*---------------------------------------------------------------------------*\
                           Class prescribedDragAve Declaration
\*---------------------------------------------------------------------------*/

class prescribedDragAve
:
    public roughnessModel
{
private:

    // Private data

        //- Residual Reynolds Number
        volScalarField dragCoeff_;
        dimensionedScalar Ct_;

public:

    //- Runtime type information
    TypeName("prescribedDragAve");


    // Constructors

        //- Construct from a dictionary and a phase pair
        prescribedDragAve
        (
            const volVectorField& U,
            const dictionary& roughnessModelDict
        );


    //- Destructor
    virtual ~prescribedDragAve();


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
