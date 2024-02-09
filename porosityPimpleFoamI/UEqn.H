// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(U)
  + 1.0/beta*fvm::div(phi, U)
  + turbulence->divDevReff(U)
 ==
    laminarTransport.nu()*fvm::laplacian(U)
  + laminarTransport.nu()/beta*(fvc::grad(beta)&fvc::grad(U)) //try times 2
  + fvm::SuSp(laminarTransport.nu()/beta*fvc::laplacian(beta),U)
  + fvOptions(beta,U)
  + pPFPorosityModel->force(U)
);

UEqn().relax();

fvOptions.constrain(UEqn());

volScalarField rAU(1.0/UEqn().A());

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
}

if(runTime.outputTime())
{ 
    volVectorField diagA("diagA", laminarTransport.nu()/beta*(fvc::grad(beta)&fvc::grad(U)));
    diagA.write();
}
