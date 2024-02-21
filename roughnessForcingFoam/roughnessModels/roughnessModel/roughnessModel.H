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
    Foam::roughnessModel

Description

SourceFiles
    roughnessModel.C
    newroughnessModel.C

\*---------------------------------------------------------------------------*/

#ifndef roughnessModel_H
#define roughnessModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "wallDist.H"
#include "fvMatricesFwd.H"
#include "fvm.H"
#include "fvc.H"
//#include "fvMatrices.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class roughnessModel Declaration
\*---------------------------------------------------------------------------*/

class roughnessModel
:
    public regIOobject
{
protected:

    // Protected data
    
        const Time& runTime_;
        const fvMesh& mesh_;
        
        const volVectorField& U_;
        
        const dictionary& roughnessModelDict_;
        
        wallDist y_;
private:

    // Private Member Functions
    
    
public:

    //- Runtime type information
    TypeName("roughnessModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            roughnessModel,
            dictionary,
            (
                const volVectorField& U,
                const dictionary& roughnessModelDict
            ),
            (U, roughnessModelDict)
        );


    // Constructors

        roughnessModel
        (
            const volVectorField& U,
            const dictionary& roughnessModelDict
        );


    //- Destructor
    virtual ~roughnessModel();


    // Selectors

        static autoPtr<roughnessModel> New
        (
            const volVectorField& U,
            const dictionary& roughnessModelDict
        );


    // Member Functions

        //- Return the roughness force used in the momentum equation
        virtual tmp<fvVectorMatrix> force(volVectorField& U) const = 0;

        //- Update roughness force for averaging
        virtual void updateForce(volVectorField& U) = 0;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
