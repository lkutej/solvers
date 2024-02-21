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

#include "myGeoUnbound.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pPFPorosityModels
{
    defineTypeNameAndDebug(myGeoUnbound, 0);
    addToRunTimeSelectionTable(pPFPorosityModel, myGeoUnbound, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

myGeoUnbound::myGeoUnbound
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

    C1_("C1", dimensionSet(0, 0, 0, 0, 0, 0, 0), pPFPorosityModelDict.lookup("C1")),

    C2_("C2", dimensionSet(0, 0, 0, 0, 0, 0, 0), pPFPorosityModelDict.lookup("C2")),

    DmMode_("DmMode", dimensionSet(0, 0, 0, 0, 0, 0, 0), pPFPorosityModelDict.lookup("DmMode")),

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
    ),

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
        dimensionedVector("fi", dimensionSet(0,1,-2,0,0,0,0), vector(0.0,0.0,0.0))
    )


{
    List<scalar> yD(pPFPorosityModelDict.lookup("y"));
    List<scalar> DmD(pPFPorosityModelDict.lookup("y"));

    if (DmMode_.value()==1.0)
    {
        DmD=pPFPorosityModelDict.lookup("Dm_ana");
        Info<<"DmMode: Dm_ana"<<endl;
        Info<<DmD<<endl;
    }
    else if (DmMode_.value()==2.0)
    {
        DmD=pPFPorosityModelDict.lookup("Dm_num");
        Info<<"DmMode: Dm_num"<<endl;
        Info<<DmD<<endl;
    }
    else if (DmMode_.value()==3.0)
    {
        DmD=pPFPorosityModelDict.lookup("sf");
        Info<<"DmMode: Dm from sf"<<endl;
        Info<<DmD<<endl;
    }

    Info<< "Interpolating Dm\n" << endl;
    scalar yL = 0.0;
    scalar yU = 0.0;
    scalar DmU = 0.0;
    scalar DmL = 0.0;
    forAll(Dm_, cellI)
    {
        //scalar y = mesh_.cellCentres()[cellI].component(1);
        forAll(yD, yI)
        {
            Dm_[cellI] = 0.0;
            //if (yD[yI] >= y)
            if (yD[yI] >= y_[cellI])
            {
                //Info<<y<<" is between "<<yD[yI-1]<<" and "<<yD[yI]<<endl;
                yL = yD[yI-1];
                yU = yD[yI];
                DmL = DmD[yI-1];
                DmU = DmD[yI];
                if (DmMode_.value()==3.0)
                {
                    Dm_[cellI] = 4.0*alpha_[cellI]/(max(DmL + (DmU-DmL)*(y_[cellI]-yL)/(yU-yL),SMALL)*constant::mathematical::pi);
                }
                else
                {
                    Dm_[cellI] = DmL + (DmU-DmL)*(y_[cellI]-yL)/(yU-yL);
                }
                break;
            }
        }
    }
/*
    if (DmMode_.value()==3.0)
    {  
        dimensionedScalar Dmin("Dmin", dimensionSet(0,1,0,0,0,0,0), SMALL); 
        Dm_ = 4.0*alpha_/(max(Dm_,Dmin)*constant::mathematical::pi);
    }
*/
    Dm_.write();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

myGeoUnbound::~myGeoUnbound()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<fvVectorMatrix> myGeoUnbound::force(volVectorField& U)
{
    dimensionedScalar Dmin("Dmin", dimensionSet(0,1,0,0,0,0,0), SMALL);
    dimensionedScalar Tmax("Tmax", dimensionSet(0,-1,0,0,0,0,0), 1000.0);
    dimensionedScalar Umin("Umin", dimensionSet(0,1,-1,0,0,0,0), SMALL);

    fi_ = -(C1_*nu_*sqr(alpha_)/(sqr((1-alpha_)*max(Dm_,Dmin))) + C2_*alpha_/((1.0-alpha_)*max(Dm_,Dmin))*mag(U))*cmptMultiply(vector(Cx_.value(), Cy_.value(), Cz_.value()),U);

    return
    (
//      fvm::SuSp(-nu*180.0*sqr(alpha_)/(sqr(d_p*(1.0-alpha_))) -nu*180.0*alpha_/(100.0*d_p*(1.0-alpha_)*nu)*mag(U),U)
//      fvm::SuSp(-(C1_*nu_*sqr(alpha_)/(sqr((1-alpha_)*max(Dm_,Dmin))) + C2_*alpha_/((1.0-alpha_)*max(Dm_,Dmin))*mag(U)),U)
      -(C1_*nu_*sqr(alpha_)/(sqr((1-alpha_)*max(Dm_,Dmin))) + C2_*alpha_/((1.0-alpha_)*max(Dm_,Dmin))*mag(U))*cmptMultiply(vector(Cx_.value(), Cy_.value(), Cz_.value()),U) + fvm::Sp(dimensionedScalar("zero", dimless/dimTime, scalar(0.0)),U)
//        -min(C1_*nu_*sqr(alpha_)/(sqr((1-alpha_)*max(Dm_,Dmin))*max(mag(U),Umin)) + C2_*alpha_/((1.0-alpha_)*max(Dm_,Dmin)),Tmax)*mag(U)*cmptMultiply(vector(Cx_.value(), Cy_.value(), Cz_.value()),U) + fvm::Sp(dimensionedScalar("zero", dimless/dimTime, scalar(0.0)),U)      
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pPFPorosityModels
} // End namespace Foam

// ************************************************************************* //
