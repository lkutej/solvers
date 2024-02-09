/*---------------------------------------------------------------------------*\
dynamicSmagorinsky - Implementation of the dynamic Smagorinsky
		     SGS model.
    
Copyright Information
    Copyright (C) 1991-2009 OpenCFD Ltd.
    Copyright (C) 2010-2014 Alberto Passalacqua 
    
License
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dynamicSmagorinsky.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace pPFLESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynamicSmagorinsky, 0);
addToRunTimeSelectionTable(pPFLESModel, dynamicSmagorinsky, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void dynamicSmagorinsky::updateSubGridScaleFields
(
    const volSymmTensorField& D
)
{
    // The SGS viscosity is bounded so that nuEff cannot become negative.
    // Values are limited here, and not in nuEff, for consistency in stored
    // data and in submodels using nuSgs().
    // No warning message is printed when this limitation is applied.
    nuSgs_ = max(cD(D)*sqr(delta())*sqrt(magSqr(D)), -nu());
    nuSgs_.correctBoundaryConditions();
}

volScalarField dynamicSmagorinsky::cD
(
    const volSymmTensorField& D
) const
{
    tmp<volSymmTensorField> LL = 
	dev(filter_(sqr(U())) - (sqr(filter_(U()))));

    const volSymmTensorField MM
    (
        sqr(delta())*(filter_(mag(D)*(D)) - 4*mag(filter_(D))*filter_(D))
    );

    // Locally averaging MMMM on cell faces
    volScalarField MMMM = fvc::average(magSqr(MM));

    MMMM.max(VSMALL);

    // Performing local average on cell faces on return
    return 0.5*fvc::average(LL && MM)/MMMM;
}


volScalarField dynamicSmagorinsky::cI
(
    const volSymmTensorField& D
) const
{
    tmp<volScalarField> KK = 
	0.5*(filter_(magSqr(U())) - magSqr(filter_(U())));

    const volScalarField mm
    (
        sqr(delta())*(4*sqr(mag(filter_(D))) - filter_(sqr(mag(D))))
    );

    // Locally averaging mmmm on cell faces
    volScalarField mmmm = fvc::average(magSqr(mm));

    mmmm.max(VSMALL);

    // Performing local average on cell faces on return
    return fvc::average(KK*mm)/mmmm;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynamicSmagorinsky::dynamicSmagorinsky
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& beta,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    pPFLESModel(typeName, U, phi, beta, transport),
    pPFGenEddyVisc(U, phi, beta, transport),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    updateSubGridScaleFields(dev(symm(fvc::grad(U))));

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dynamicSmagorinsky::correct
(
    const tmp<volTensorField>& gradU
)
{
    pPFLESModel::correct(gradU);

    const volSymmTensorField D(dev(symm(gradU)));

    k_ = cI(D)*sqr(delta())*magSqr(D);
    bound(k_,  kMin_);

    updateSubGridScaleFields(D);
}

bool dynamicSmagorinsky::read()
{
    if (pPFGenEddyVisc::read())
    {
        filter_.read(coeffDict());

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
