/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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

#include "deLemoskEp.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
//#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace pPFRASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(deLemoskEp, 0);
addToRunTimeSelectionTable(pPFRASModel, deLemoskEp, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> deLemoskEp::fMu() const
{
    return exp(-3.4/sqr(scalar(1) + sqr(k_)/(nu()*epsilon_)/50.0));
}


tmp<volScalarField> deLemoskEp::f2() const
{
    return
        scalar(1)
      - 0.3*exp(-min(sqr(sqr(k_)/(nu()*epsilon_)), scalar(50.0)));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

deLemoskEp::deLemoskEp
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& beta,
    transportModel& transport,
    const word& pPFTurbulenceModelName,
    const word& modelName
)
:
    pPFRASModel(modelName, U, phi, beta, transport, pPFTurbulenceModelName),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    CEps1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps1",
            coeffDict_,
            1.44
        )
    ),
    CEps2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            coeffDict_,
            1.9
        )
    ),
    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK",
            coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.33
        )
    ),
    ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            0.28
        )
    ),
 
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

    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
   
    yr_(mesh_)
    
{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);
    printCoeffs();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<volSymmTensorField> deLemoskEp::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> deLemoskEp::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> deLemoskEp::divDevReff(volVectorField& U) const
{
    return
    (
        -1.0/beta_*
        (
            fvm::laplacian(nut_*beta_, U)
          + fvc::div(nut_*U*fvc::grad(beta_))
          + fvc::div(nut_*dev2(T(fvc::grad(beta_*U))))
        )
    );
}

tmp<fvVectorMatrix> deLemoskEp::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}

bool deLemoskEp::read()
{
    if (pPFRASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        CEps1_.readIfPresent(coeffDict());
        CEps2_.readIfPresent(coeffDict());
        sigmaK_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        ck_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void deLemoskEp::correct()
{
    pPFRASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    volScalarField G("pPFRASModel::G", 1.0/beta_*nut_*2*magSqr(symm(fvc::grad(beta_*U_)))); 
    //laut de Lemos -R:grad(beta*U). Da -beta R definiert ist mit Boussinesq, hier 1/beta_.
    
    /* Porosity source term */
    dimensionedScalar dp("dp", dimLength, 0.01);
    volScalarField K("K", pow(dp,2.0)*pow(beta_,3.0)/(180*pow(max(1.0-beta_,SMALL),2.0)));
    volScalarField B("B", ck_*beta_*k_*mag(beta_*U_)/sqrt(K)); 
    Info<<"Source term, ck = "<<ck_<<endl;

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(beta_, epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff()*beta_, epsilon_)
      - fvc::laplacian(DepsilonEff()*epsilon_, beta_)
     ==
        CEps1_*G/k_*epsilon_
      - fvm::Sp(f2()*CEps2_/k_*epsilon_*beta_, epsilon_)
      + fvm::Sp(CEps2_*ck_*beta_*mag(beta_*U_)/sqrt(K), epsilon_)
    );

    epsEqn().relax();
    //epsEqn().boundaryManipulate(epsilon_.boundaryField());

    dimensionedScalar nu1 = nu()->internalField()[0];
    #include "BoundaryConditionEp.H"

    solve(epsEqn);
    bound(epsilon_, epsilonMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(beta_, k_) //beta
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff()*beta_, k_) //beta
      - fvc::laplacian(DkEff()*k_, beta_) //added
     ==
        G
      - fvm::Sp(beta_*epsilon_/k_, k_) //beta
      + ck_*beta_*k_*mag(beta_*U_)/sqrt(K) //Source term added
     );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);
        
    // Re-calculate viscosity
    nut_ = Cmu_*fMu()*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();

    Info<<"k, min: "<<min(k_).value()<<" max: "<<max(k_).value()<<" average: "<<k_.weightedAverage(mesh_.V()).value()<<endl;
    Info<<"epsilon, min: "<<min(epsilon_).value()<<" max: "<<max(epsilon_).value()<<" average: "<<epsilon_.weightedAverage(mesh_.V()).value()<<endl;

    if(runTime_.outputTime())
    {  
        B.write();
        K.write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pPFRASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
