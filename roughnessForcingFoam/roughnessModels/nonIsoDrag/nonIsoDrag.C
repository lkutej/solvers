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

#include "nonIsoDrag.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace roughnessModels
{
    defineTypeNameAndDebug(nonIsoDrag, 0);
    addToRunTimeSelectionTable(roughnessModel, nonIsoDrag, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nonIsoDrag::nonIsoDrag
(
    const volVectorField& U,
    const dictionary& roughnessModelDict
)
:
    roughnessModel(U, roughnessModelDict),
    
    dragCoeff_
    (
        IOobject
        (
            "dragCoeff",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ), 

    Ct_("Ct", dimless, roughnessModelDict.lookup("Ct")),

    fi_
    (
        IOobject
        (
            "fi",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("fi", dimensionSet(0, 1, -2, 0, 0, 0, 0), vector::zero)
    )
{
    List<scalar> yD(roughnessModelDict.lookup("y"));
    List<scalar> dragCoeffD(roughnessModelDict.lookup("dragCoeff"));

    Info<< "Interpolating dragCoeff\n" << endl;
    scalar yL = 0.0;
    scalar yU = 0.0;
    scalar dragCoeffU = 0.0;
    scalar dragCoeffL = 0.0;
    forAll(dragCoeff_, cellI)
    {
        //scalar y = mesh_.cellCentres()[cellI].component(1);
        forAll(yD, yI)
        {
            dragCoeff_[cellI] = 0.0;
            //if (yD[yI] >= y)
            if (yD[yI] >= y_[cellI])
            {
                //Info<<y<<" is between "<<yD[yI-1]<<" and "<<yD[yI]<<endl;
                yL = yD[yI-1];
                yU = yD[yI];
                dragCoeffL = dragCoeffD[yI-1];
                dragCoeffU = dragCoeffD[yI];
                dragCoeff_[cellI] = dragCoeffL + (dragCoeffU-dragCoeffL)*(y_[cellI]-yL)/(yU-yL);
                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

nonIsoDrag::~nonIsoDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<fvVectorMatrix> nonIsoDrag::force(volVectorField& U) const
{
    dimensionedScalar dimFix("dimFix", dimless/dimLength, 1.0);
    
    return
    (
        fvm::SuSp(-dragCoeff_ * Ct_ * dimFix * mag(U), U) 
        /*
        -0.5*Cd_*1.0/k_*max((1.0-y_/k_),0.0)*cmptMultiply(vector(1.0, 0.0, 0.0),sign(U.component(0))*cmptMultiply(U,U)) 
        + fvm::Sp(dimensionedScalar("zero", dimless/dimTime, scalar(0.0)),U)
        */
    );
}

void nonIsoDrag::updateForce(volVectorField& U)
{
    dimensionedScalar dimFix("dimFix", dimless/dimLength, 1.0);
    fi_ = -dragCoeff_ * Ct_ * dimFix * mag(U) * U;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace roughnessModels
} // End namespace Foam

// ************************************************************************* //
