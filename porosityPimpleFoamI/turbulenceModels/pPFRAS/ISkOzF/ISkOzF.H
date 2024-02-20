/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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
    Foam::incompressible::RASModels::ISkOzF

Description
    Zeta-effe turbulence model for incompressible flows.
    
    Reference article: 

    "A robust near-wall elliptic-relaxation eddy-viscosity turbulence model for CFD"
    M. Popovac et al. (2004), International Journal of Heat and Fluid Flow 25 (2004)

    The default model coefficients correspond to the following:
    \verbatim
        ISkOzFCoeffs
        {
             Cmu_         0.22,
             CEps2_       1.90,
             C1_          0.40,
             C2_          0.65,
             sigmaK_      1.00,
             sigmaEps_    1.30,
             sigmaZ_      1.20,
             CTau_        6.00,
             CL_          0.36,
             CEta_        85.0
        }
    \endverbatim

SourceFiles
    ISkOzF.C

\*---------------------------------------------------------------------------*/

#ifndef ISkOzF_H
#define ISkOzF_H

#include "RASModel.H"
#include "wallDistReflection.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class ISkOzF Declaration
\*---------------------------------------------------------------------------*/

class ISkOzF
:
    public RASModel
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar CEps2_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar sigmaCDv_;
            dimensionedScalar sigmaCDt_;
            dimensionedScalar sigmaZ_;
            dimensionedScalar CTau_;
            dimensionedScalar CL_;
            dimensionedScalar CEta_;
            dimensionedScalar a_;
            dimensionedScalar TSwitch_;
            dimensionedScalar LSwitch_;
            dimensionedScalar Csas_;
            dimensionedScalar CT2_;
            dimensionedScalar Clim_;
            dimensionedScalar fwSwitch_;
            dimensionedScalar fEqnPbykSwitch_;
            dimensionedScalar DavidsonSwitch_;
            dimensionedScalar LTaylorSwitch_;
            dimensionedScalar SdiffSwitch_;
            dimensionedScalar CDtboundSwitch_;
            dimensionedScalar Ceps1zetaSwitch_;

        // Fields

            volScalarField k_;
            volScalarField kr_;
            volScalarField epsilon_;
            volScalarField epsilonr_;
            volScalarField omega_;
            volScalarField zeta_;
            volScalarField f_;
            volScalarField nut_;
            wallDistReflection yr_;
            volScalarField delta_;
            volScalarField T1_;
            volScalarField T2_;

       // averaging
            IOdictionary VLESAveragingProperties_;
            dimensionedScalar Dt_;
            volVectorField UAvg_;

       //- Lower limit of zeta
 	    dimensionedScalar mTSmall_;
            dimensionedScalar zetaMin_;

       //- Lower limit of f
            dimensionedScalar fMin_;

       // Private member functions

            tmp<volScalarField> Tau() const;
            tmp<volScalarField> L() const;
            void calculateDelta();
            void writeAveragingProperties() const;
	// Wall distance reflection
	    

public:

    //- Runtime type information
    TypeName("ISkOzF");

    // Constructors

        //- Construct from components
        ISkOzF
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~ISkOzF()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the lower allowable limit for zeta (default: SMALL)
        const dimensionedScalar& zetaMin() const
        {
            return zetaMin_;
        }

        //- Allow zetaMin to be changed
        dimensionedScalar& zetaMin()
        {
            return zetaMin_;
        }

        //- Return the lower allowable limit for f (default: SMALL)
        const dimensionedScalar& fMin() const
        {
            return fMin_;
        }

        //- Allow fMin to be changed
        dimensionedScalar& fMin()
        {
            return fMin_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_/sigmaK_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + nu())
            );
        }

        //- Return the effective diffusivity for zeta
        tmp<volScalarField> DzetaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DzetaEff", nut_/sigmaZ_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the v2-k ratio (zeta)
        virtual tmp<volScalarField> zeta() const
        {
            return zeta_;
        }

        virtual tmp<volScalarField> f() const
        {
            return f_;
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
