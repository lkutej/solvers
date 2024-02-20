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

#include "VLES.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
//#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(VLES, 0);
addToRunTimeSelectionTable(RASModel, VLES, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> VLES::Tau() const
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
    Info<<"T.Org, "; 
    return max
           (
               min
               (
                   k_/(epsilon_+epsilonMin_),
                   a_/((sqrt(6.0)*Cmu_*mag(symm(fvc::grad(U_))))*zeta_+mTSmall_)    //BK: Durbin Original 1/(sqrt(6) zeta Cmu sqrt(Sij Sij)
               ),
               CTau_*sqrt(nu()/(epsilon_+epsilonMin_))
           ); 
*/
    
    //Info<<"T.Org without realizability constraint, "; 

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

/*
    Info<<"Tau = turbulent time scale, "; 
    return k_/epsilon_;
*/
/*    
    return max(k_/(epsilon_+epsilonMin_), sqrt(nu()/(epsilon_+epsilonMin_)));
*/
}

tmp<volScalarField> VLES::L() const
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
    Info<<"L.org"<<endl;
    return CL_*max
               (
                   min                                                              //BK: in Samules Version auskommentiert
                   (                                                                //BK: in Samules Version auskommentiert
                       pow(k_,1.5)/(epsilon_+epsilonMin_),
                       sqrt(k_)/((sqrt(6.0)*Cmu_*mag(symm(fvc::grad(U_))))*zeta_)   //BK: in Samules Version auskommentiert
                   ),                                                               //BK: in Samules Version auskommentiert
                   CEta_*pow( (pow(nu(),3)/(epsilon_+epsilonMin_)),0.25)
               );
*/

    //Lt_ = CL_ * pow(k_,1.5)/(epsilon_+epsilonMin_);
    //Lk_ = CL_ * CEta_ * pow((pow(nu(),3)/(epsilon_+epsilonMin_)),0.25);
    
    volScalarField Li = CL_ * pow(k_,1.5)/epsilon_;
    volScalarField Lk = CL_ * LSwitch_*CEta_ * pow(pow(nu(),3)/epsilon_,0.25);
    volScalarField LInd("LInd", pos(Li-Lk));

    if (runTime_.outputTime())
    {
        LInd.write();
    }

    //Info<<"L.org without realizability constraint"<<endl;
    return CL_*max
               (
                   pow(k_,1.5)/epsilon_,
                   LSwitch_*CEta_*pow(pow(nu(),3)/epsilon_,0.25)
               );

/*	       
    return CL_*max(pow(k_,1.5)/(epsilon_+epsilonMin_), pow((pow(nu(),3)/(epsilon_+epsilonMin_)),0.25));
*/
}

void VLES::calculateDelta()
{
    /*
    // ~~~~~~~~~~~~~~~~~~~~~~~ max(dx,dy,dz) ~~~~~~~~~~~~~~~~~~~~~~~ //
  
    Info<<"delta=max(dx,dy,dz)"<<endl;

    tmp<volScalarField> hmax
    (
        new volScalarField
        (
            IOobject
            (
                "hmax",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zrero", dimLength, 0.0)
        )
    );

    const cellList& cells = mesh_.cells();
    const vectorField& cellC = mesh_.cellCentres();
    const vectorField& faceC = mesh_.faceCentres();
    const vectorField faceN(mesh_.faceAreas()/mag(mesh_.faceAreas()));

    forAll(cells, cellI)
    {
        scalar deltaMaxTmp = 0.0;
        const labelList& cFaces = cells[cellI];
        const point& cc = cellC[cellI];

        forAll(cFaces, cFaceI)
        {
            label faceI = cFaces[cFaceI];
            const point& fc = faceC[faceI];
            const vector& n = faceN[faceI];

            scalar tmp = magSqr(n*(n & (fc - cc)));
            if (tmp > deltaMaxTmp)
            {
                deltaMaxTmp = tmp;
            }
        }

        hmax()[cellI] = 2.0 * sqrt(deltaMaxTmp);
    }

    delta_.internalField() = hmax();
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    */
    
    
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~ mean  ~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    /* 
    Info<<"delta = 1/3*(dx+dy+dz)"<<endl;
    
    tmp<volScalarField> hsum
    (
        new volScalarField
        (
            IOobject
            (
                "hsum",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zrero", dimLength, 0.0)
        )
    );

    const cellList& cells = mesh_.cells();
    const vectorField& cellC = mesh_.cellCentres();
    const vectorField& faceC = mesh_.faceCentres();
    const vectorField faceN(mesh_.faceAreas()/mag(mesh_.faceAreas()));

    forAll(cells, cellI)
    {
        const labelList& cFaces = cells[cellI];
        const point& cc = cellC[cellI];

        int ii=0;
        forAll(cFaces, cFaceI)
        {
            label faceI = cFaces[cFaceI];
            const point& fc = faceC[faceI];
            const vector& n = faceN[faceI];

            scalar tmp = magSqr(n*(n & (fc - cc)));
            hsum()[cellI] += sqrt(tmp);
            ii++;
        }
	Info<<ii<<endl;
    }

    delta_.internalField() = 1.0/3.0 * hsum();
    */
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    
    
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~ cube root ~~~~~~~~~~~~~~~~~~~~~~~~~ //
    
    Info<<"delta = mesh.V^(1/3)"<<endl;
    delta_.internalField() = pow(mesh_.V(), 1.0/3.0);
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    
}
	
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

VLES::VLES
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

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

    ksgs_
    (
        IOobject
        (
            "ksgs",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("ksgs", dimensionSet(0, 2, -2, 0, 0, 0, 0), 1.0e-10)
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
    
    fk_
    (
        IOobject
        (
            "fk",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_, 0.1
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
 
    yr_(mesh_),
    
    mTSmall_("mTSmall", dimless/dimTime, SMALL),
    //kMin_("kMin", sqr(dimVelocity), 1.0e-10),
    //epsilonMin_("epsilonMin", kMin_.dimensions()/dimTime, SMALL),
    zetaMin_("zetaMin", dimless, SMALL),
    fMin_("fMin", dimless/dimTime, SMALL),
    TscMin_("TscMin", dimTime, SMALL),
    LscMin_("LscMin", dimLength, SMALL)
{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);
    bound(f_, fMin_);
    bound(zeta_, zetaMin_);

    nut_ = Cmu_*zeta_*sqr(fk_)*k_*Tau();
    //nut_.correctBoundaryConditions();
   
    calculateDelta();    

    printCoeffs();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<volSymmTensorField> VLES::R() const
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


tmp<volSymmTensorField> VLES::devReff() const
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


tmp<fvVectorMatrix> VLES::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}

tmp<fvVectorMatrix> VLES::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev2(T(fvc::grad(U))))
    );
}

bool VLES::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        CEps2_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmaK_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        sigmaZ_.readIfPresent(coeffDict());
        CTau_.readIfPresent(coeffDict());
        CL_.readIfPresent(coeffDict());
        CEta_.readIfPresent(coeffDict());
        a_.readIfPresent(coeffDict());
        TSwitch_.readIfPresent(coeffDict());
        LSwitch_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void VLES::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    volScalarField G("RASModel::G", nut_ * 2 * magSqr(symm(fvc::grad(U_))));
    //volScalarField CEps1_ = 1.4*(1.0+(0.012/zeta_));
    volScalarField CEps1_ = 1.4*(1.0+0.045/sqrt(zeta_)); 
    
    /* 
    //~~~~~~~~~~~~~ update length and time scale ~~~~~~~~~~~~~
    Tt_ = k_/epsilon_;
    Tk_ = CTau_*sqrt(nu()/epsilon_);
    Lt_ = CL_ * pow(k_,1.5)/epsilon_;
    Lk_ = CL_ * CEta_*pow(pow(nu(),3)/epsilon_,0.25);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */    

    volScalarField T_ = Tau();// + TscMin_;
    volScalarField L_ = L();// + LscMin_;
    
    // Update epsilon (and possibly G) at the wall
    // epsilon_.boundaryField().updateCoeffs();
     
    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        CEps1_*G/T_
      - fvm::Sp(CEps2_/T_, epsilon_)
/*
        CEps1_*G*epsilon_/k_
      - fvm::Sp(CEps2_*epsilon_/k_, epsilon_)
*/
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
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::Sp(epsilon_/k_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // f equation
    tmp<fvScalarMatrix> fEqn
    (
        - fvm::laplacian(f_)
     ==
        - fvm::Sp(1.0/sqr(L_),f_)
        - (C1_+C2_*G/epsilon_)*(zeta_ - 2.0/3.0)/(sqr(L_)*T_)
//        - (C1_+C2_*G/epsilon_)*(zeta_ - 2.0/3.0)/(sqr(L_)*k_/epsilon_)
    );
    
    fEqn().relax();
    solve(fEqn);
    bound(f_, fMin_);

    
    // Calculate f from the transformed fEqn
    volScalarField fTilda = f_ - 2.0*nu()*zeta_/sqr(yr_);
    
    
    // Zeta equation
    tmp<fvScalarMatrix> zetaEqn
    (
        fvm::ddt(zeta_)
      + fvm::div(phi_, zeta_)
      - fvm::laplacian(DzetaEff(), zeta_)
     ==
//      f_
        fTilda
      - fvm::Sp(G/(k_), zeta_)
    );

    zetaEqn().relax();
    solve(zetaEqn);
    bound(zeta_, zetaMin_);
    //zeta_ = min(zeta_, 2.0);

    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    if (mesh_.changing())
    {
	calculateDelta();
    };

    fk_ = max(min(pow(delta_/(pow(k_, 1.5)/epsilon_), 2.0/3.0), 1.0), 1.0e-5);
    volScalarField Fr("Fr", sqr(fk_));
      
    ksgs_ = fk_*k_;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    // Re-calculate viscosity
    nut_ = Cmu_*zeta_*k_*T_*Fr;
    //nut_.correctBoundaryConditions();

    //omega calc for VLESkOSST
    volScalarField omega("omega", epsilon_/(0.09*k_));
    
    if(runTime_.outputTime())
    {
        Fr.write();
        omega.write();
        volScalarField LvK("LvK", 0.41*sqrt(2.0*magSqr(symm(fvc::grad(U_))))/(mag(fvc::laplacian(U_))));
        LvK.write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
