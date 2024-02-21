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

#include "prescribedDragAveTS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace roughnessModels
{
    defineTypeNameAndDebug(prescribedDragAveTS, 0);
    addToRunTimeSelectionTable(roughnessModel, prescribedDragAveTS, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

prescribedDragAveTS::prescribedDragAveTS
(
    const volVectorField& U,
    const dictionary& roughnessModelDict
)
:
    roughnessModel(U, roughnessModelDict),
    
    dragCoeff1_
    (
        IOobject
        (
            "dragCoeff1",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ), 

    dragCoeff2_
    (
        IOobject
        (
            "dragCoeff2",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    C1_("C1", dimless, roughnessModelDict.lookup("C1")),
 
    C2_("C2", dimless, roughnessModelDict.lookup("C2")),

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
    List<scalar> dragCoeffD1(roughnessModelDict.lookup("dragCoeff1"));
    List<scalar> dragCoeffD2(roughnessModelDict.lookup("dragCoeff2"));

    Info<< "Interpolating dragCoeff\n" << endl;
    scalar yL = 0.0;
    scalar yU = 0.0;
    vector yC;
    scalar dragCoeffU = 0.0;
    scalar dragCoeffL = 0.0;
    scalar trapz1 = 0.0;
    scalar trapz2 = 0.0;

    const faceList& ff = mesh_.faces();
    const pointField& pp = mesh_.points();

    forAll(y_, cellI)
    {
        const cell& cc = mesh_.cells()[cellI];
        labelList pLabels(cc.labels(ff));
        //Info<<"pLabels: "<<pLabels<<endl;
        pointField pLocal(pLabels.size(), vector::zero);

        forAll (pLabels, pointI)
        {
            pLocal[pointI] = pp[pLabels[pointI]];
            //Info<<"pLocal[pointI]: "<<pLocal[pointI]<<endl;
        }

        yL = min(pLocal & vector(0,1,0));
        yU = max(pLocal & vector(0,1,0));

        trapz1 = 0.0;
        trapz2 = 0.0;

        forAll(yD, iI)
        {
            if ((yD[iI] > yL) && (yD[iI] < yU))
            {
                //Info<<"Point "<<yD[iI]<<" is in this cell"<<endl;
                if (yD[iI-1] <= yL) // If true, first point in cell
                {
                    // Interpolate value at yL and add segement from yL to current point to trapz
                    dragCoeffL = dragCoeffD1[iI-1] + (dragCoeffD1[iI]-dragCoeffD1[iI-1])*(yL-yD[iI-1])/(yD[iI]-yD[iI-1]);
                    trapz1 = trapz1 + (yD[iI]-yL)*0.5*(dragCoeffL+dragCoeffD1[iI]);
                    dragCoeffL = dragCoeffD2[iI-1] + (dragCoeffD2[iI]-dragCoeffD2[iI-1])*(yL-yD[iI-1])/(yD[iI]-yD[iI-1]);
                    trapz2 = trapz2 + (yD[iI]-yL)*0.5*(dragCoeffL+dragCoeffD2[iI]);
                }
                else
                {
                    // Add segment from previous point to point to trapz
                    trapz1 = trapz1 + (yD[iI]-yD[iI-1])*0.5*(dragCoeffD1[iI-1]+dragCoeffD1[iI]);
                    trapz2 = trapz2 + (yD[iI]-yD[iI-1])*0.5*(dragCoeffD2[iI-1]+dragCoeffD2[iI]);
                }
                if (yD[iI+1] > yU) // If true, last point in cell
                {
                    // Interpolate value at yU and add segment from current point to yU to trapz
                    dragCoeffU = dragCoeffD1[iI] + (dragCoeffD1[iI+1]-dragCoeffD1[iI])*(yU-yD[iI])/(yD[iI+1]-yD[iI]);
                    trapz1 = trapz1 + (yU-yD[iI])*0.5*(dragCoeffD1[iI]+dragCoeffU);
                    dragCoeffU = dragCoeffD2[iI] + (dragCoeffD2[iI+1]-dragCoeffD2[iI])*(yU-yD[iI])/(yD[iI+1]-yD[iI]);
                    trapz2 = trapz2 + (yU-yD[iI])*0.5*(dragCoeffD2[iI]+dragCoeffU);
                }
            }
            else if ((yD[iI] > yU) && (yD[iI-1] < yL)) // No point in this cell
            {
                // Interpolate value at yL and yU and add segment to trapz
                dragCoeffL = dragCoeffD1[iI-1] + (dragCoeffD1[iI]-dragCoeffD1[iI-1])*(yL-yD[iI-1])/(yD[iI]-yD[iI-1]);
                dragCoeffU = dragCoeffD1[iI-1] + (dragCoeffD1[iI]-dragCoeffD1[iI-1])*(yU-yD[iI-1])/(yD[iI]-yD[iI-1]);
                trapz1 = trapz1 + (yU-yL)*0.5*(dragCoeffL+dragCoeffU);
                dragCoeffL = dragCoeffD2[iI-1] + (dragCoeffD2[iI]-dragCoeffD2[iI-1])*(yL-yD[iI-1])/(yD[iI]-yD[iI-1]);
                dragCoeffU = dragCoeffD2[iI-1] + (dragCoeffD2[iI]-dragCoeffD2[iI-1])*(yU-yD[iI-1])/(yD[iI]-yD[iI-1]);
                trapz2 = trapz2 + (yU-yL)*0.5*(dragCoeffL+dragCoeffU);
            }
        }
        dragCoeff1_[cellI] = trapz1/(yU-yL);
        dragCoeff2_[cellI] = trapz2/(yU-yL);
    }
    dragCoeff1_.write();
    dragCoeff2_.write();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

prescribedDragAveTS::~prescribedDragAveTS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<fvVectorMatrix> prescribedDragAveTS::force(volVectorField& U) const
{
    dimensionedScalar dimFix1("dimFix1", dimless/dimTime, 1.0);
    dimensionedScalar dimFix2("dimFix2", dimless/dimLength, 1.0);
    
    return
    (
        fvm::SuSp(-(C1_*dragCoeff1_*dimFix1+C2_*dragCoeff2_*dimFix2*mag(U)), U) 
        /*
        -0.5*Cd_*1.0/k_*max((1.0-y_/k_),0.0)*cmptMultiply(vector(1.0, 0.0, 0.0),sign(U.component(0))*cmptMultiply(U,U)) 
        + fvm::Sp(dimensionedScalar("zero", dimless/dimTime, scalar(0.0)),U)
        */
    );
}

void prescribedDragAveTS::updateForce(volVectorField& U)
{
    dimensionedScalar dimFix1("dimFix1", dimless/dimTime, 1.0);
    dimensionedScalar dimFix2("dimFix2", dimless/dimLength, 1.0);
    fi_ = -(C1_*dragCoeff1_*dimFix1+C2_*dragCoeff2_*dimFix2*mag(U))*U;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace roughnessModels
} // End namespace Foam

// ************************************************************************* //
