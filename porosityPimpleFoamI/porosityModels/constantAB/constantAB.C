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

#include "constantAB.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pPFPorosityModels
{
    defineTypeNameAndDebug(constantAB, 0);
    addToRunTimeSelectionTable(pPFPorosityModel, constantAB, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constantAB::constantAB
(
    const volVectorField& U,
    const volScalarField& alpha,
    const dictionary& pPFPorosityModelDict
)
:
    pPFPorosityModel(U, alpha, pPFPorosityModelDict),
    A_("A", dimensionSet(0, 0, 0, 0, 0, 0, 0), pPFPorosityModelDict.lookup("A")),
    B_("B", dimensionSet(0, 0, 0, 0, 0, 0, 0), pPFPorosityModelDict.lookup("B"))

{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constantAB::~constantAB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<fvVectorMatrix> constantAB::force(volVectorField& U) const
{
   dimensionedScalar dimFixA("dimFix", dimless/dimTime, 1.0);
   dimensionedScalar dimFixB("dimFix", dimless/dimLength, 1.0);
   // volScalarField dragCoeff("dragCoeff", Cd_*pos(k_-y_)*dimFix);
    return
    (
        fvm::SuSp(-pos(alpha_-SMALL)*(A_*dimFixA + B_*dimFixB*mag(U)),U)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pPFPorosityModels
} // End namespace Foam

// ************************************************************************* //
