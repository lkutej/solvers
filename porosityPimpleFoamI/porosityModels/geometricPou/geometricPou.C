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

#include "geometricPou.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pPFPorosityModels
{
    defineTypeNameAndDebug(geometricPou, 0);
    addToRunTimeSelectionTable(pPFPorosityModel, geometricPou, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

geometricPou::geometricPou
(
    const volVectorField& U,
    const volScalarField& alpha,
    const dictionary& pPFPorosityModelDict
)
:
    pPFPorosityModel(U, alpha, pPFPorosityModelDict),

    Cx_("Cx", dimensionSet(0, 0, 0, 0, 0, 0, 0), pPFPorosityModelDict.lookup("Cx")),

    Cy_("Cy", dimensionSet(0, 0, 0, 0, 0, 0, 0), pPFPorosityModelDict.lookup("Cy")),

    Cz_("Cz", dimensionSet(0, 0, 0, 0, 0, 0, 0), pPFPorosityModelDict.lookup("Cz")),

    nu_("nu", dimensionSet(0, 2, -1, 0, 0, 0, 0), pPFPorosityModelDict.lookup("nu")),

    CD_("CD", dimensionSet(0, 0, 0, 0, 0, 0, 0), pPFPorosityModelDict.lookup("CD")),

    CF_("CF", dimensionSet(0, 0, 0, 0, 0, 0, 0), pPFPorosityModelDict.lookup("CF")),

    s_
    (
        IOobject
        (
            "s",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("s", dimensionSet(0,-1,0,0,0,0,0), 0.0)
    ),

    sf_
    (
        IOobject
        (
            "sf",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("sf", dimensionSet(0,-1,0,0,0,0,0), 0.0)
    )

{
    List<scalar> yD(pPFPorosityModelDict.lookup("y"));
    List<scalar> sD(pPFPorosityModelDict.lookup("y"));
    List<scalar> sfD(pPFPorosityModelDict.lookup("y"));

    sD=pPFPorosityModelDict.lookup("s");
    sfD=pPFPorosityModelDict.lookup("sf");

    Info<< "Interpolating s\n" << endl;
    scalar yL = 0.0;
    scalar yU = 0.0;
    scalar vU = 0.0;
    scalar vL = 0.0;
    forAll(s_, cellI)
    {
        //scalar y = mesh_.cellCentres()[cellI].component(1);
        forAll(yD, yI)
        {
            s_[cellI] = 0.0;
            //if (yD[yI] >= y)
            if (yD[yI] >= y_[cellI])
            {
                //Info<<y<<" is between "<<yD[yI-1]<<" and "<<yD[yI]<<endl;
                yL = yD[yI-1];
                yU = yD[yI];
                vL = sD[yI-1];
                vU = sD[yI];
                s_[cellI] = vL + (vU-vL)*(y_[cellI]-yL)/(yU-yL);
                break;
            }
        }
    }
    s_.write();

    Info<< "Interpolating sf\n" << endl;
    yL = 0.0;
    yU = 0.0;
    vU = 0.0;
    vL = 0.0;
    forAll(sf_, cellI)
    {
        //scalar y = mesh_.cellCentres()[cellI].component(1);
        forAll(yD, yI)
        {
            sf_[cellI] = 0.0;
            //if (yD[yI] >= y)
            if (yD[yI] >= y_[cellI])
            {
                //Info<<y<<" is between "<<yD[yI-1]<<" and "<<yD[yI]<<endl;
                yL = yD[yI-1];
                yU = yD[yI];
                vL = sfD[yI-1];
                vU = sfD[yI];
                s_[cellI] = vL + (vU-vL)*(y_[cellI]-yL)/(yU-yL);
                break;
            }
        }
    }
    sf_.write();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

geometricPou::~geometricPou()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<fvVectorMatrix> geometricPou::force(volVectorField& U) const
{
    if(runTime_.outputTime())
    {
        volVectorField fiGeo("fiGeo", -(CD_*nu_*sqr(s_)/pow(max(1.0-alpha_,SMALL),3.0) + CF_*sf_/2.0*mag(U))*cmptMultiply(vector(Cx_.value(), Cy_.value(), Cz_.value()),U));
        fiGeo.write();
    }

    return
    (
        //fvm::SuSp(-(nu_*2.0*C1_*alpha_/(constant::mathematical::pi*max((1.0-alpha_)*sqr(Dm_),sqr(Dmin))) + 2.0*C2_*alpha_/(constant::mathematical::pi*max((1.0-alpha_)*Dm_,Dmin))*mag(U)),U)
        -(CD_*nu_*sqr(s_)/pow(max(1.0-alpha_,SMALL),3.0) + CF_*sf_/2.0*mag(U))*cmptMultiply(vector(Cx_.value(), Cy_.value(), Cz_.value()),U) + fvm::Sp(dimensionedScalar("zero", dimless/dimTime, scalar(0.0)),U)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pPFPorosityModels
} // End namespace Foam

// ************************************************************************* //
