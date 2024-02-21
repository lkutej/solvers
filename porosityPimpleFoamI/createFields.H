Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//#include "createPhi.H"

// *************** BK: reading alpha and calculating phi ***************

Info<< "Reading field alpha\n" << endl;
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField beta("beta", 1.0-alpha);

Info<< "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(beta*U) & mesh.Sf()
);

// *********************************************************************


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::pPFTurbulenceModel> turbulence
(
    incompressible::pPFTurbulenceModel::New(U, phi, beta, laminarTransport)
);

// ********* BK: reading porosityModelDict and creating model *********

IOdictionary porosityModelDict
(
    IOobject
    (
        "porosityModelDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

autoPtr<pPFPorosityModel> pPFPorosityModel = pPFPorosityModel::New
(
    U,
    alpha,
    porosityModelDict
);

// *********************************************************************

volSymmTensorField Rmod
(
    IOobject
    (   
        "Rmod",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    turbulence->R()
);

volVectorField Us
(
    IOobject
    (   
        "Us",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U*beta
);
