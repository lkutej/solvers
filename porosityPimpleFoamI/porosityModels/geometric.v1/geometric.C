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

#include "geometric.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pPFPorosityModels
{
    defineTypeNameAndDebug(geometric, 0);
    addToRunTimeSelectionTable(pPFPorosityModel, geometric, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

geometric::geometric
(
    const volVectorField& U,
    const volScalarField& alpha,
    const dictionary& pPFPorosityModelDict
)
:
    pPFPorosityModel(U, alpha, pPFPorosityModelDict),

    nu_("nu", dimensionSet(0, 2, -1, 0, 0, 0, 0), pPFPorosityModelDict.lookup("nu")),

    C1_("C1", dimensionSet(0, 0, 0, 0, 0, 0, 0), pPFPorosityModelDict.lookup("C1")),

    C2_("C2", dimensionSet(0, 0, 0, 0, 0, 0, 0), pPFPorosityModelDict.lookup("C2")),

    N_("N", dimensionSet(0, 0, 0, 0, 0, 0, 0), pPFPorosityModelDict.lookup("N")),

    Dm_
    (
        IOobject
        (
            "Dm",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Dm", dimensionSet(0,1,0,0,0,0,0), 0.0)
    )

{
    dimensionedScalar As("As", dimensionSet(0,2,0,0,0,0,0), 32.0);
    Dm_ = sqrt(4.0*alpha_*As/(constant::mathematical::pi*N_));
    Dm_.write();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

geometric::~geometric()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<fvVectorMatrix> geometric::force(volVectorField& U) const
{
    dimensionedScalar Dmin("Dmin", dimensionSet(0,1,0,0,0,0,0), SMALL);

    if(runTime_.outputTime())
    {
        volVectorField fiGeo("fiGeo", -(nu_*2.0*C1_*alpha_/(constant::mathematical::pi*max((1.0-alpha_)*sqr(Dm_),sqr(Dmin))) + 2.0*C2_*alpha_/(constant::mathematical::pi*max((1.0-alpha_)*Dm_,Dmin))*mag(U))*U);
        fiGeo.write();
    }

    return
    (
//      fvm::SuSp(-nu*180.0*sqr(alpha_)/(sqr(d_p*(1.0-alpha_))) -nu*180.0*alpha_/(100.0*d_p*(1.0-alpha_)*nu)*mag(U),U)
        fvm::SuSp(-(nu_*2.0*C1_*alpha_/(constant::mathematical::pi*max((1.0-alpha_)*sqr(Dm_),sqr(Dmin))) + 2.0*C2_*alpha_/(constant::mathematical::pi*max((1.0-alpha_)*Dm_,Dmin))*mag(U)),U)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pPFPorosityModels
} // End namespace Foam

// ************************************************************************* //
