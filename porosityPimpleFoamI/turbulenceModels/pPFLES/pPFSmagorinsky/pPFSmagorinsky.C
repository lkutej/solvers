/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "pPFSmagorinsky.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace pPFLESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pPFSmagorinsky, 0);
addToRunTimeSelectionTable(pPFLESModel, pPFSmagorinsky, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void pPFSmagorinsky::updateSubGridScaleFields(const volTensorField& gradU)
{
    nuSgs_ = ck_*delta()*sqrt(k(gradU));
    nuSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pPFSmagorinsky::pPFSmagorinsky
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& beta,
    transportModel& transport,
    const word& pPFTurbulenceModelName,
    const word& modelName
)
:
    pPFLESModel(modelName, U, phi, beta, transport, pPFTurbulenceModelName),
    pPFGenEddyVisc(U, phi, beta, transport),

    ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            0.094
        )
    )
{
    updateSubGridScaleFields(fvc::grad(U));

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pPFSmagorinsky::correct(const tmp<volTensorField>& gradU)
{
    pPFGenEddyVisc::correct(gradU);
    updateSubGridScaleFields(gradU());
}


bool pPFSmagorinsky::read()
{
    if (pPFGenEddyVisc::read())
    {
        ck_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pPFLESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
