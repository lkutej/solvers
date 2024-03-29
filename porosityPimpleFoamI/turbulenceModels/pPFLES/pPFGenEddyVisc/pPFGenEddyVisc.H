/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Class
    Foam::incompressible::LESModels::pPFGenEddyVisc

Group
    grpIcoLESTurbulence

Description
    General base class for all incompressible models that can be implemented
    as an eddy viscosity, i.e. algebraic and one-equation models.

    Contains fields for k (SGS turbulent kinetic energy), gamma
    (modelled viscosity) and epsilon (SGS dissipation).

SourceFiles
    pPFGenEddyVisc.C

\*---------------------------------------------------------------------------*/

#ifndef pPFGenEddyVisc_H
#define pPFGenEddyVisc_H

#include "pPFLESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace pPFLESModels
{

/*---------------------------------------------------------------------------*\
                         Class pPFGenEddyVisc Declaration
\*---------------------------------------------------------------------------*/

class pPFGenEddyVisc
:
    virtual public pPFLESModel
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        pPFGenEddyVisc(const pPFGenEddyVisc&);
        pPFGenEddyVisc& operator=(const pPFGenEddyVisc&);


protected:

    // Model coefficients

        dimensionedScalar ce_;

    // Fields

        volScalarField nuSgs_;


public:

    //- Partial Runtime type information
    static const word typeName;

    // Constructors

        //- Construct from components
        pPFGenEddyVisc
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& beta,
            transportModel& transport,
            const word& pPFTurbulenceModelName = pPFTurbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~pPFGenEddyVisc()
    {}


    // Member Functions

        //- Return sub-grid disipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (   
                        "epsilon",
                        runTime_.timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    ce_*k()*sqrt(k())/delta()
                )
            );
        }

        //- Return the SGS viscosity.
        virtual tmp<volScalarField> nuSgs() const
        {
            return nuSgs_;
        }

        //- Return the sub-grid stress tensor.
        virtual tmp<volSymmTensorField> B() const;

        //- Return the effective sub-grid turbulence stress tensor
        //  including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pPFLESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
