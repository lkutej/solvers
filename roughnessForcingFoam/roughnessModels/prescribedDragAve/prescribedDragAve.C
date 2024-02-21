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

#include "prescribedDragAve.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace roughnessModels
{
    defineTypeNameAndDebug(prescribedDragAve, 0);
    addToRunTimeSelectionTable(roughnessModel, prescribedDragAve, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

prescribedDragAve::prescribedDragAve
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

    Ct_("Ct", dimless, roughnessModelDict.lookup("Ct"))
{
    List<scalar> yD(roughnessModelDict.lookup("y"));
    List<scalar> dragCoeffD(roughnessModelDict.lookup("dragCoeff"));

    // This is only correct for hex-cells!
    Info<< "Interpolating dragCoeff\n" << endl;
    scalar yL = 0.0;
    scalar yU = 0.0;
    vector yC;
    scalar dragCoeffL = 0.0;
    scalar dragCoeffU = 0.0;
    scalar trapz = 0.0;

    // new test:------------------------------------------------------------------------------
    const faceList& ff = mesh_.faces();
    const pointField& pp = mesh_.points();
    // new test end --------------------------------------------------------------------------

    forAll(dragCoeff_, cellI)
    {
        /*
        // This works, if cells are only in y direction
        yL = yU;
        yU = yL + 2.0 * (y_[cellI]-yU);
        Info<<"Cell "<<cellI<<" goes from y="<<yL<<" to y="<<yU<<endl;
        */

        // new test: --------------------------------------------------------------------------
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

        // new test end -----------------------------------------------------------------------

        trapz = 0.0;
        forAll(yD, iI)
        {
            if ((yD[iI] > yL) && (yD[iI] < yU))
            {
                //Info<<"Point "<<yD[iI]<<" is in this cell"<<endl;
                if (yD[iI-1] <= yL) // If true, first point in cell
                {
                    // Interpolate value at yL and add segement from yL to current point to trapz
                    dragCoeffL = dragCoeffD[iI-1] + (dragCoeffD[iI]-dragCoeffD[iI-1])*(yL-yD[iI-1])/(yD[iI]-yD[iI-1]);
                    trapz = trapz + (yD[iI]-yL)*0.5*(dragCoeffL+dragCoeffD[iI]);
                    //Info<<"Point "<<yD[iI]<<" is first point in this cell"<<endl;
                    //Info<<"trapz is "<<trapz<<endl;
                }
                else
                {
                    // Add segment from previous point to point to trapz
                    trapz = trapz + (yD[iI]-yD[iI-1])*0.5*(dragCoeffD[iI-1]+dragCoeffD[iI]);
                    //Info<<"trapz is "<<trapz<<endl;
                }
                if (yD[iI+1] > yU) // If true, last point in cell
                {
                    // Interpolate value at yU and add segment from current point to yU to trapz
                    dragCoeffU = dragCoeffD[iI] + (dragCoeffD[iI+1]-dragCoeffD[iI])*(yU-yD[iI])/(yD[iI+1]-yD[iI]);
                    trapz = trapz + (yU-yD[iI])*0.5*(dragCoeffD[iI]+dragCoeffU);
                    // Maybe break here?
                    //Info<<"Point "<<yD[iI]<<" is last point in this cell"<<endl;
                    //Info<<"trapz is "<<trapz<<endl;
                }
            }
            else if ((yD[iI] > yU) && (yD[iI-1] < yL)) // No point in this cell
            {
                // Interpolate value at yL and yU and add segment to trapz
                dragCoeffL = dragCoeffD[iI-1] + (dragCoeffD[iI]-dragCoeffD[iI-1])*(yL-yD[iI-1])/(yD[iI]-yD[iI-1]);
                dragCoeffU = dragCoeffD[iI-1] + (dragCoeffD[iI]-dragCoeffD[iI-1])*(yU-yD[iI-1])/(yD[iI]-yD[iI-1]);
                trapz = trapz + (yU-yL)*0.5*(dragCoeffL+dragCoeffU);
            }
        }
        //Info<<"trapz is "<<trapz<<endl;
        dragCoeff_[cellI] = trapz/(yU-yL);
        //Info<<"dragCoeff_[cellI] is "<<dragCoeff_[cellI]<<endl;
    }
    dragCoeff_.write();


/*
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
*/
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

prescribedDragAve::~prescribedDragAve()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<fvVectorMatrix> prescribedDragAve::force(volVectorField& U) const
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace roughnessModels
} // End namespace Foam

// ************************************************************************* //
