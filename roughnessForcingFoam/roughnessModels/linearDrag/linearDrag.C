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

#include "linearDrag.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace roughnessModels
{
    defineTypeNameAndDebug(linearDrag, 0);
    addToRunTimeSelectionTable(roughnessModel, linearDrag, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearDrag::linearDrag
(
    const volVectorField& U,
    const dictionary& roughnessModelDict
)
:
    roughnessModel(U, roughnessModelDict),
    Cd_("Cd", dimless, roughnessModelDict.lookup("Cd")),
    k_("k", dimLength, roughnessModelDict.lookup("k")),
    r0_("r0", dimless, roughnessModelDict.lookup("r0")),
    r1_("r1", dimless, roughnessModelDict.lookup("r1")),
    alpha_
    (
        IOobject
        (
            "alpha",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

linearDrag::~linearDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<fvVectorMatrix> linearDrag::force(volVectorField& U) const
{
    volScalarField yWall = y_;
    dimensionedScalar dimFix("dimFix", dimless/dimLength, 1.0);
    volScalarField dragCoeff("dragCoeff", Cd_ * min((r0_+(r1_-r0_)*y_/k_),1.0) * pos(k_-y_) * dimFix);

    if(runTime_.outputTime())
    {
        yWall.write();
        dragCoeff.write();
    }
    return
    (
        fvm::SuSp(-dragCoeff*mag(U),U)
        /*
        -0.5*Cd_*1.0/k_*max((1.0-y_/k_),0.0)*cmptMultiply(vector(1.0, 0.0, 0.0),sign(U.component(0))*cmptMultiply(U,U)) 
        + fvm::Sp(dimensionedScalar("zero", dimless/dimTime, scalar(0.0)),U)
        */
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace roughnessModels
} // End namespace Foam

// ************************************************************************* //
