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

#include "kOzFp.H"
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

defineTypeNameAndDebug(kOzFp, 0);
addToRunTimeSelectionTable(pPFRASModel, kOzFp, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> kOzFp::Tau() const
{
/*    
    volScalarField T_lb("T_lb", CTau_*sqrt(nu()/(epsilon_+epsilonMin_)));
    volScalarField T_ub("T_ub", a_/((sqrt(6.0)*Cmu_*mag(symm(fvc::grad(U_))))*zeta_+mTSmall_));
    volScalarField T_nb("T_nb",  k_/(epsilon_+epsilonMin_));
    
    dimensionedScalar TSmall("TSmall", dimTime, 1e-15);
    
    volScalarField v_min = min(T_nb, T_ub);
    volScalarField I_min = 1.0*pos(T_nb - T_ub - TSmall); // = 1 wenn gebounded
    volScalarField I_max = 2.0*pos(T_lb - v_min - TSmall); // = 1 wenn gebounded
    volScalarField TInd("TInd", I_min + I_max);
    
    if (runTime_.outputTime())
    {
        TInd.write();
        T_lb.write();
        T_ub.write();
        T_nb.write();
    }
*/ 
/*  
    return max(
                  min(
                           k_/(epsilon_+epsilonMin_),
                                 ( a_/((sqrt(6.0)*Cmu_*mag(symm(fvc::grad(U_))))*zeta_+mTSmall_))
                       ),
                   CTau_*sqrt(nu()/(epsilon_+epsilonMin_))
               );
*/

    volScalarField Ti = k_/epsilon_;
    volScalarField Tk = TSwitch_*CTau_ * sqrt(nu()/epsilon_);
    volScalarField TInd("TInd", pos(Ti-Tk));

    if (runTime_.outputTime())
    {
        TInd.write();
    }
    
    return max
           (
               k_/epsilon_,
               TSwitch_*CTau_*sqrt(nu()/epsilon_)
           );
}

tmp<volScalarField> kOzFp::L() const
{
/*
    volScalarField L_lb("L_lb", CEta_*pow( (pow(nu(),3)/(epsilon_+epsilonMin_)),0.25));
    volScalarField L_ub("L_ub",sqrt(k_)/((sqrt(6.0)*Cmu_*mag(symm(fvc::grad(U_))))*zeta_));
    volScalarField L_nb("L_nb",pow(k_,1.5)/(epsilon_+epsilonMin_));
    
    dimensionedScalar LSmall("LSmall", dimLength, 1e-15);
    
    volScalarField v_min = min(L_nb, L_ub);
    volScalarField I_min = 1.0*pos(L_nb - L_ub - LSmall); // = 1 wenn gebounded
    volScalarField I_max = 2.0*pos(L_lb - v_min - LSmall); // = 1 wenn gebounded
    volScalarField LInd("LInd", I_min + I_max);
    
    if (runTime_.outputTime())
    {
        LInd.write();
        L_lb.write();
        L_ub.write();
        L_nb.write();
    }
*/
/*
    return CL_*max(
                      min(
                              pow(k_,1.5)/(epsilon_+epsilonMin_),
                                sqrt(k_)/((sqrt(6.0)*Cmu_*mag(symm(fvc::grad(U_))))*zeta_)
                          ),
                      CEta_*pow( (pow(nu(),3)/(epsilon_+epsilonMin_)),0.25)
                  );
*/

    volScalarField Li = CL_ * pow(k_,1.5)/epsilon_;
    volScalarField Lk = CL_ * LSwitch_*CEta_ * pow(pow(nu(),3)/epsilon_,0.25);
    volScalarField LInd("LInd", pos(Li-Lk));

    if (runTime_.outputTime())
    {
        LInd.write();
    }

    if (LTaylorSwitch_.value() == 1)
    {
        Info<<"USING 13*L_TAYLOR"<<endl;
        return 13.0*sqrt(10.0*nu()*k_/epsilon_);
    }
    else
    {
        return CL_*max
                 (
                     pow(k_,1.5)/epsilon_,
                     LSwitch_*CEta_*pow(pow(nu(),3)/epsilon_,0.25)
                 );
    }
/*
    return CL_*max
               (
                   pow(k_,1.5)/epsilon_,
                   LSwitch_*CEta_*pow(pow(nu(),3)/epsilon_,0.25)
               );
*/
/*
    Info<<"USING 13*L_TAYLOR"<<endl;
    return 13.0*sqrt(10.0*nu()*k_/epsilon_);
*/
}

void kOzFp::calculateDelta()
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~ cube root ~~~~~~~~~~~~~~~~~~~~~~~~~ //

    Info<<"delta = mesh.V^(1/3)"<<endl;
    delta_.internalField() = pow(mesh_.V(), 1.0/3.0);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
}

void kOzFp::writeAveragingProperties() const
{
    IOdictionary propsDict
    (
        IOobject
        (
            "VLESAveragingProperties",
            mesh_.time().timeName(),
            "uniform",
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    propsDict.add("Dt", Dt_);
    propsDict.regIOobject::write();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kOzFp::kOzFp
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
            0.22
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
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            0.4
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            0.65
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
            1.3
        )
    ),
    sigmaCDv_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaCDv",
            coeffDict_,
            1.4
        )
    ),
    sigmaCDt_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaCDt",
            coeffDict_,
            1.5
        )
    ),
    sigmaZ_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaZ",
            coeffDict_,
            1.2
        )
    ),
    CTau_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTau",
            coeffDict_,
            6.0
        )
    ),
    CL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CL",
            coeffDict_,
            0.36
        )
    ),
    CEta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CEta",
            coeffDict_,
            85
        )
    ),
    a_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a",
            coeffDict_,
            0.6
        )
    ),
    TSwitch_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "TSwitch",
            coeffDict_,
            0.0
        )
    ),
    LSwitch_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "LSwitch",
            coeffDict_,
            0.0
        )
    ),
    Csas_
    (   
        dimensioned<scalar>::lookupOrAddToDict
        (   
            "Csas",
            coeffDict_,
            3.0
        )
    ),
    CT2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CT2",
            coeffDict_,
            20.0
        )
    ),
    Clim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clim",
            coeffDict_,
            0.0
        )
    ),
    fwSwitch_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fwSwitch",
            coeffDict_,
            0.0
        )
    ), 
    fEqnPbykSwitch_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fEqnPbykSwitch",
            coeffDict_,
            0.0
        )
    ),
    DavidsonSwitch_
    (   
        dimensioned<scalar>::lookupOrAddToDict
        (   
            "DavidsonSwitch",
            coeffDict_,
            0.0
        )
    ),
    LTaylorSwitch_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "LTaylorSwitch",
            coeffDict_,
            0.0
        )
    ),
    SdiffSwitch_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "SdiffSwitch",
            coeffDict_,
            1.0
        )
    ),
    CDtboundSwitch_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDtboundSwitch",
            coeffDict_,
            0.0
        )
    ),
    Ceps1zetaSwitch_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps1zetaSwitch",
            coeffDict_,
            1.0
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

    co_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "co",
            coeffDict_,
            1.0
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

    kr_
    (
        IOobject
        (
            "kr",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("kr", dimensionSet(0, 2, -2, 0, 0, 0, 0), 1.0e-10)
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

    epsilonr_
    (
        IOobject
        (
            "epsilonr",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("epsilonr", dimensionSet(0, 2, -3, 0, 0, 0, 0), 0.0)
    ),

    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    zeta_
    (
        IOobject
        (
            "zeta",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    f_
    (
        IOobject
        (
            "f",
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
   
    yr_(mesh_),
    
    delta_
    (
        IOobject
        (
            "delta",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("delta", dimLength, 1.0e-10)
    ),

    T1_
    (
        IOobject
        (
            "T1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("T1", dimless/sqr(dimTime), 0.0)
    ),

    T2_
    (
        IOobject
        (
            "T2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("T2", dimless/sqr(dimTime), 0.0)
    ),

    VLESAveragingProperties_
    (
        IOobject
        (
            "VLESAveragingProperties",
            mesh_.time().timeName(),
            "uniform",
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),

    Dt_("Dt", dimless, VLESAveragingProperties_.lookup("Dt")),

    UAvg_
    (
        IOobject
        (
            "UAvg",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("UAvg", dimensionSet(0,1,-1,0,0,0,0), vector::zero)
    ), 

    mTSmall_
    (
	"mTSmall",	
	dimensionSet(0, 0, -1, 0, 0, 0, 0),
	1e-10
    ),
	
    zetaMin_
    (
        "zetaMin",
         dimless,
         SMALL
    ),
    fMin_
    (
        "fMin",
        dimless/dimTime,
        SMALL
    )
{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);
    bound(zeta_, zetaMin_);
    calculateDelta();
    printCoeffs();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<volSymmTensorField> kOzFp::R() const
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


tmp<volSymmTensorField> kOzFp::devReff() const
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


tmp<fvVectorMatrix> kOzFp::divDevReff(volVectorField& U) const
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

tmp<fvVectorMatrix> kOzFp::divDevRhoReff
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

bool kOzFp::read()
{
    if (pPFRASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        CEps2_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmaK_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        sigmaCDv_.readIfPresent(coeffDict());
        sigmaCDt_.readIfPresent(coeffDict());
        sigmaZ_.readIfPresent(coeffDict());
        CTau_.readIfPresent(coeffDict());
        CL_.readIfPresent(coeffDict());
        CEta_.readIfPresent(coeffDict());
        a_.readIfPresent(coeffDict());
        TSwitch_.readIfPresent(coeffDict());
        LSwitch_.readIfPresent(coeffDict());
        Csas_.readIfPresent(coeffDict());
        CT2_.readIfPresent(coeffDict());
        Clim_.readIfPresent(coeffDict());
        fwSwitch_.readIfPresent(coeffDict());
        fEqnPbykSwitch_.readIfPresent(coeffDict());
        DavidsonSwitch_.readIfPresent(coeffDict());
        LTaylorSwitch_.readIfPresent(coeffDict());
        SdiffSwitch_.readIfPresent(coeffDict());
        CDtboundSwitch_.readIfPresent(coeffDict());
        Ceps1zetaSwitch_.readIfPresent(coeffDict());
        ck_.readIfPresent(coeffDict());
        co_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void kOzFp::correct()
{
    pPFRASModel::correct();

    if (!turbulence_)
    {
        return;
    }


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ online averaging  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    dimensionedScalar dt = runTime_.deltaTValue();
    Dt_ += dt;

    Info<<"Dt_ is: "<<Dt_<<endl;

    dimensionedScalar alpha = (Dt_ - dt)/Dt_;
    dimensionedScalar beta  = dt/Dt_;

    //RAvg_ += sqr(UAvg_);
    
    UAvg_ = alpha*UAvg_ + beta*U_;

    //RAvg_ = alpha*RAvg_ + beta*sqr(U_) - sqr(UAvg_);

    kr_ = 0.5 * magSqr(U_-UAvg_);

    volTensorField graduPrime = fvc::grad(U_-UAvg_);
    volTensorField epsilonrTens_ = 2*nu()*(graduPrime.T() & graduPrime);
    epsilonr_ = 0.5*tr(epsilonrTens_);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    volScalarField G("pPFRASModel::G", 1.0/beta_*nut_*2*magSqr(symm(fvc::grad(beta_*U_)))); 
    //laut de Lemos -R:grad(beta*U). Da -beta R definiert ist mit Boussinesq, hier 1/beta_.
    
    volScalarField CEps1_ = 0.0*zeta_;
    if (Ceps1zetaSwitch_.value() == 1.0)
    {   
        CEps1_ = 1.4*(1.0+(0.012/(zeta_+zetaMin_)));
        Info<<"Ceps1 with zeta"<<endl;
    }
    else if (Ceps1zetaSwitch_.value() == 0.0)
    {
        CEps1_ = 1.4252+0.0*zeta_;
        Info<<"Ceps1 w/o zeta"<<endl;
    }
   
    /* Porosity source term */
    dimensionedScalar dp("dp", dimLength, 0.01);
    volScalarField K("K", pow(dp,2.0)*pow(beta_,3.0)/(180*pow(max(1.0-beta_,SMALL),2.0)));
    volScalarField F("F", beta_*dp/(100.0*max(1.0-beta_,SMALL))); //eigentlich noch 1/nu() fuer Ftilda
    //volScalarField B("B", ck_*beta_*k_*mag(beta_*U_)/sqrt(K)); //deLemos
    //volScalarField B("B", -2.0*beta_*(nu()/K+F/K*mag(U_))); //Lee
    //volScalarField kInf("kInf", pos(1.0-beta_-SMALL)*3.7*(1.0-beta_)*pow(beta_,1.5)*pow(mag(U_),2.0));
    //volScalarField epsInf("epsInf", pos(1.0-beta_-SMALL)*39.0*sqr(beta_)*pow(1.0-beta_,2.5)*1.0/dp*pow(mag(U_),3.0));
    Info<<"Source term deLemos"<<endl;
    if(runTime_.outputTime())
    {
        //epsInf.write();
        //kInf.write();
        //dimSwitch.write();
        //K.write();
        //B.write();
    }

    dimensionedScalar dimFix("dimFix", dimensionSet(0,-1,0,0,0,0,0), 1.0);
    volScalarField  dimSwitch("dimSwitch", pos(1.0-beta_-SMALL)*dimFix);

    volScalarField T_ = Tau();
    volScalarField L_ = L();

    volScalarField CDv("CDv", (2.0/k_*nu()/sigmaCDv_*(fvc::grad(k_)&fvc::grad(omega_)))/omega_);
    volScalarField CDt("CDt", (2.0/k_*nut_/sigmaCDt_*(fvc::grad(k_)&fvc::grad(omega_)))/omega_);
    volScalarField CD("CD", CDv + CDt);
    if (CDtboundSwitch_.value() == 1.0)
    {
        CD=CDv+max(CDt,dimensionedScalar("0.0", dimless/dimTime, 0.0));
        Info<<"Bounded turbulent cross diffusion"<<endl;
    }

    Info<<"derivative beta terms in omegaEqn disabled"<<endl;
 
    // omega equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(beta_, omega_) //beta
      + fvm::div(phi_, omega_)
      - fvm::laplacian(DepsilonEff()*beta_, omega_) //beta
      - fvc::laplacian(DepsilonEff()*omega_, beta_) //added
      ==
        (CEps1_-1.0)*G/k_*omega_
      - fvm::SuSp((CEps2_-1.0)*omega_*beta_, omega_) //beta
      + fvm::SuSp(co_*beta_*mag(beta_*U_)/sqrt(K), omega_) //deLemosSource
//    + fvm::SuSp((CEps2_-1.0)*co_*beta_*mag(U_)*F/K, omega_) //deLemosSourceFK
//    + fvm::SuSp((CEps2_-1.0)*co_*mag(U_)*dimSwitch, omega_) //deLemosSourceCali
//    + fvm::SuSp(-2.0*co_*beta_*(nu()/K+F/K*mag(U_)), omega_) //LeeSource
//    + beta_/k_*(3.0*CEps2_*sqr(epsInf)/max(kInf,kMin_)-epsInf*omega_) //NakayamaSource
//    + fvm::SuSp(fvc::div(phi_, beta_), omega_) //Transformation term
//    - fvc::laplacian(DkEff(), beta_)*omega_ //Transformation term
      + fvm::Sp(beta_*CD, omega_) // Cross diffusion //beta
      + fvm::SuSp(SdiffSwitch_*(fvc::laplacian(DepsilonEff(), k_) - fvc::laplacian(DkEff(), k_))/k_*beta_, omega_) //Zero, if sigmaEps=sigmaK //beta
    );
    omegaEqn().relax();
    dimensionedScalar nu1 = nu()->internalField()[0];
    #include "BoundaryConditionOmega.H"
    solve(omegaEqn);
    bound(omega_, omegaMin_);
    //omega_.max(pol);

    epsilon_ = omega_ * k_;
    bound(epsilon_, epsilonMin_);

/*
    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        CEps1_*G/T_
      - fvm::Sp(CEps2_/T_, epsilon_)
      + Csas_*Psas*k_
    );

    epsEqn().relax();
    //epsEqn().boundaryManipulate(epsilon_.boundaryField());

    dimensionedScalar nu1 = nu()->internalField()[0];
    #include "BoundaryConditionEp.H"

    solve(epsEqn);
    bound(epsilon_, epsilonMin_);
*/

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
      + ck_*beta_*mag(beta_*U_)/sqrt(K)*k_ //deLemosSource
//    + ck_*beta_*mag(U_)*F/K*k_ //deLemosSourceFK
//    + ck_*mag(U_)*dimSwitch*k_ //deLemosSourceFK
//    - 2.0*beta_*(nu()/K+ck_*F/K*mag(U_))*k_ //LeeSource
//    + beta_*epsInf //NakayamaSource
     );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);

    if (fwSwitch_.value() == 0.0)
    {
        Info<<"Using fTilda"<<endl;
        // f equation
        tmp<fvScalarMatrix> fEqn
        (
            - fvm::laplacian(beta_,f_) //beta
         ==
            - fvm::Sp(beta_/sqr(L_),f_) //beta
//          - (C1_+C2_*G/epsilon_)*(zeta_ - 2.0/3.0)/(sqr(L_)*T_)
            - beta_/(sqr(L_)*T_)*C1_*(zeta_-2.0/3.0) //beta
            - 1.0/(sqr(L_)*T_)*C2_*G/epsilon_*(zeta_-2.0/3.0)
            + fEqnPbykSwitch_ * 1.0/sqr(L_)*2.0/3.0*G/k_
        );
        
        fEqn().relax();
        solve(fEqn);
        bound(f_, fMin_);

        //Calculate f from the transformed fEqn
        volScalarField fTilda = f_ - 2.0*nu()*zeta_/sqr(yr_);

        // zeta equation
        tmp<fvScalarMatrix> zetaEqn
        (
            fvm::ddt(beta_, zeta_) //beta
          + fvm::div(phi_, zeta_)
          - fvm::laplacian(DzetaEff()*beta_, zeta_) //beta
          - fvc::laplacian(DzetaEff()*zeta_, beta_) //added
         ==
            beta_*fTilda //beta
          - fvm::Sp(G/k_, zeta_)
        );

        zetaEqn().relax();
        solve(zetaEqn);
        bound(zeta_, zetaMin_);
        zeta_ = min(zeta_, 2.0);
    }
    else if (fwSwitch_.value() == 1.0)
    {
        Info<<"Using f"<<endl;
        /*  
        T_ = Tau();
        L_ = L();
        */
        // f equation
        tmp<fvScalarMatrix> fEqn
        (
            - fvm::laplacian(beta_, f_) //beta
         ==
            - fvm::Sp(beta_/sqr(L_),f_) //beta
//          - (C1_+(C2_*G/(epsilon_)))*((zeta_ - 2.0/3.0))/(sqr(L_)*T_)
            - beta_/(sqr(L_)*T_)*C1_*(zeta_-2.0/3.0) //beta
            - 1.0/(sqr(L_)*T_)*C2_*G/epsilon_*(zeta_-2.0/3.0)
            + fEqnPbykSwitch_ * 1.0/sqr(L_)*2.0/3.0*G/k_
        );
        
        fEqn().relax();
        #include "BoundaryConditionf.H"
        solve(fEqn);

        // zeta equation
        tmp<fvScalarMatrix> zetaEqn
        (
            fvm::ddt(beta_,zeta_) //beta
          + fvm::div(phi_, zeta_)
          - fvm::laplacian(DzetaEff()*beta_, zeta_) //beta
          - fvc::laplacian(DzetaEff()*zeta_, beta_) //added
         ==
//            min(f_, (C1_+(C2_*G/(epsilon_)))*((zeta_ - 2.0/3.0))/T_)
            beta_*f_ //beta
          - fvm::Sp(G/k_, zeta_)
        );

        zetaEqn().relax();
        solve(zetaEqn);
        bound(zeta_, zetaMin_);
        zeta_ = min(zeta_, 2.0);
    }
        
    // Re-calculate viscosity
    if (DavidsonSwitch_.value() == 0.0)
    {
        //nut_ = pos(beta_-1.0)*Cmu_*zeta_*k_*T_;
        //Info<<"nut disabled in porous region!"<<endl;        
        nut_ = Cmu_*zeta_*k_*T_;
    }
    else if (DavidsonSwitch_.value() == 1.0)
    {
        nut_ = min(Cmu_*zeta_*k_*T_, 0.09*sqr(k_)/epsilon_);
    }
    else if (DavidsonSwitch_.value() == 2.0)
    {
        nut_ = 0.09*sqr(k_)/epsilon_;
        Info<<"kE-mode"<<endl;
    }
    nut_.correctBoundaryConditions();

    Info<<"k, min: "<<min(k_).value()<<" max: "<<max(k_).value()<<" average: "<<k_.weightedAverage(mesh_.V()).value()<<endl;
    Info<<"epsilon, min: "<<min(epsilon_).value()<<" max: "<<max(epsilon_).value()<<" average: "<<epsilon_.weightedAverage(mesh_.V()).value()<<endl;

    if(runTime_.outputTime())
    {  
        writeAveragingProperties();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pPFRASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
