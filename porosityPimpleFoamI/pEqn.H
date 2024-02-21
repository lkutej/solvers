surfaceScalarField rAUf("rAUf", fvc::interpolate(beta*rAU));

volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn().H();
/*
if (pimple.nCorrPISO() <= 1)
{
    UEqn.clear();
}
*/
surfaceScalarField phiHbyA
(
    "phiHbyA",
    (fvc::interpolate(beta*HbyA) & mesh.Sf())
  + rAUf*fvc::ddtCorr(beta, U, phi) //()?
);

fvOptions.makeRelative(fvc::interpolate(beta), phiHbyA); //()?

adjustPhi(phiHbyA, U, p); //()?

// Update the fixedFluxPressure BCs to ensure flux consistency
setSnGrad<fixedFluxPressureFvPatchScalarField> //()?
(
    p.boundaryField(),
    (
        phiHbyA.boundaryField()
      - fvOptions.relative(mesh.Sf().boundaryField() & U.boundaryField())
       *beta.boundaryField() //()?
    )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}
/*
if(runTime.outputTime())
{ 
    volVectorField resBeta("resBeta", fvc::reconstruct(phi)-beta*U);
    resBeta.write();
    volVectorField recPhi("recPhi", fvc::reconstruct(phi));
    recPhi.write();
    volScalarField contiV("contiV", fvc::div(beta*U));
    contiV.write();
    volScalarField contiP("contiP", fvc::div(phi));
    contiP.write();
}
*/
#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

Rmod = turbulence->R(); //BK: Update Rmod
Us = U*beta; //BK: Update superficial velocity
