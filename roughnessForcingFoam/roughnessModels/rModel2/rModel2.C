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

#include "rModel2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace roughnessModels
{
    defineTypeNameAndDebug(rModel2, 0);
    addToRunTimeSelectionTable(roughnessModel, rModel2, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rModel2::rModel2
(
    const volVectorField& U,
    const dictionary& roughnessModelDict
)
:
    roughnessModel(U, roughnessModelDict),
    Cd_("Cd", dimless, roughnessModelDict.lookup("Cd")),
    k_("k", dimLength, roughnessModelDict.lookup("k")),
    AfbyAt_("AfbyAt", dimless, roughnessModelDict.lookup("AfbyAt")),
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

rModel2::~rModel2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<fvVectorMatrix> rModel2::force(volVectorField& U) const
{
    volScalarField yWall = y_;
    volScalarField dragCoeff("dragCoeff", 0.5*Cd_/(1-alpha_)*1.0/k_*AfbyAt_*pos(k_-y_));


    if(runTime_.outputTime())
    {
        yWall.write();
        dragCoeff.write();
    }
    //dimensionedScalar dimFix("dimFix", dimensionSet(0,-1,0,0,0,0,0), 1.0);
    //volScalarField dragCoeff(y_);
    Info<<k_<<endl;
    return
    (
        fvm::SuSp(-dragCoeff*mag(U),U)
        /*
        -0.5*Cd_*1.0/k_*max((1.0-y_/k_),0.0)*cmptMultiply(vector(1.0, 0.0, 0.0),sign(U.component(0))*cmptMultiply(U,U)) 
        + fvm::Sp(dimensionedScalar("zero", dimless/dimTime, scalar(0.0)),U)
        */
        //fvm::SuSp(mag(U),U
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace roughnessModels
} // End namespace Foam

// ************************************************************************* //
