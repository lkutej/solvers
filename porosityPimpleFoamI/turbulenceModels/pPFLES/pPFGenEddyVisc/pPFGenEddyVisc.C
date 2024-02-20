/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "pPFGenEddyVisc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace pPFLESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameWithName(pPFGenEddyVisc, "pPFGenEddyVisc");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pPFGenEddyVisc::pPFGenEddyVisc
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

    ce_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ce",
            coeffDict_,
            1.048
        )
    ),

    nuSgs_
    (
        IOobject
        (
            "nuSgs",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
//    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> pPFGenEddyVisc::B() const
{
    return ((2.0/3.0)*I)*k() - nuSgs_*twoSymm(fvc::grad(U()));
}


tmp<volSymmTensorField> pPFGenEddyVisc::devReff() const
{
    return -nuEff()*dev(twoSymm(fvc::grad(U())));
}


tmp<fvVectorMatrix> pPFGenEddyVisc::divDevReff(volVectorField& U) const
{
    return
    (
      -1.0/beta_*
      (
          fvm::laplacian(nuSgs_*beta_, U)
          + fvc::div(nuSgs_*U*fvc::grad(beta_))
          + fvc::div(nuSgs_*dev2(T(fvc::grad(beta_*U))))
      )
    );
}


tmp<fvVectorMatrix> pPFGenEddyVisc::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}


void pPFGenEddyVisc::correct(const tmp<volTensorField>& gradU)
{
    pPFLESModel::correct(gradU);
}


bool pPFGenEddyVisc::read()
{
    if (pPFLESModel::read())
    {
        ce_.readIfPresent(coeffDict());

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
