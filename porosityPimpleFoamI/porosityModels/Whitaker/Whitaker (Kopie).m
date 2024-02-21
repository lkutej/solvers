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

\*---------------------------------------------------------------------------*/

#include "Whitaker.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pPFPorosityModels
{
    defineTypeNameAndDebug(Whitaker, 0);
    addToRunTimeSelectionTable(pPFPorosityModel, Whitaker, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Whitaker::Whitaker
(
    const volVectorField& U,
    const volScalarField& alpha,
    const dictionary& pPFPorosityModelDict
)
:
    pPFPorosityModel(U, alpha, pPFPorosityModelDict),
    nu("nu", dimensionSet(0, 2, -1, 0, 0, 0, 0), pPFPorosityModelDict.lookup("nu")),
    d_p("d_p", dimensionSet(0, 1, 0, 0, 0, 0, 0), pPFPorosityModelDict.lookup("d_p"))

{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Whitaker::~Whitaker()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<fvVectorMatrix> Whitaker::force(volVectorField& U) const
{
    volScalarField yWall = y_;
   //dimensionedScalar dimFix("dimFix", dimless/dimLength, 1.0);
   // volScalarField dragCoeff("dragCoeff", Cd_*pos(k_-y_)*dimFix);
 

    if(runTime_.outputTime())
    {
        yWall.write();
     //  dragCoeff.write();

  
    }


    return
    (

	  fvm::SuSp(-nu*((180.0*sqr(alpha_))/(sqr(d_p)*(1.0-alpha_)))*(1.0+(((1.0-alpha_)*d_p)/(100.0*alpha_*nu)*mag(U))),U)  - fvm::laplacian((1.0-alpha_),U)
      //fvm::SuSp(-0.5*mag(U),U)
		
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pPFPorosityModels
} // End namespace Foam

// ************************************************************************* //
