Info<< "Reading TransportProperties\n" << endl;

 #include "readTransportProperties.H"

Info<< "Reading field enthalpy\n" << endl;
    volScalarField enthalpy
    (
        IOobject
        (
            "enthalpy",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
Info<< "Reading field T\n" << endl;

     volScalarField T
     (
    IOobject
     (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
      ),
         Foam::pos(cpS*Ts-enthalpy)*(enthalpy/cpS)
      	+Foam::pos(enthalpy-(cpS*Ts + hs))*(Tl+(enthalpy-(cpS*Ts + hs))/cpL)
      	+Foam::pos(enthalpy-cpS*Ts)*Foam::pos((cpS*Ts + hs)-enthalpy)*(Ts+ (enthalpy-cpS*Ts)*(Tl-Ts)/hs)
      );

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
        max(    min((T-Ts)/(Tl-Ts),scalar(1))   ,scalar(0))
    );
Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
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

#include "createPhi.H"


//Info<< "Creating turbulence model\n" << endl;
//autoPtr<incompressible::turbulenceModel> turbulence
//(
//    incompressible::turbulenceModel::New(U, phi, laminarTransport)
//);

// Kinematic density for buoyancy force
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh
    ),
  //  1.0 - Foam::max(beta*(T - Tl), 0.0)
      1.0 - Foam::max(beta*(T - Ts), 0.0)

);

    // D'arcy-type source term field
    volScalarField DC
    (
        IOobject
        (
            "DC",
            runTime.timeName(),
            mesh
        ),
        DCl*Foam::pow(1.0-alpha,2)/(Foam::pow(alpha,3)+DCs)
    );

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh
        ),
        alpha*rhoL+(1.0-alpha)*rhoS
    );

 //   surfaceScalarField rhophi
 //    (
 //   IOobject
  //   (
   //     "rhophi",
   //     runTime.timeName(),
   //     mesh,
  //      IOobject::READ_IF_PRESENT,
  //      IOobject::AUTO_WRITE
  //    ),
  //  linearInterpolate(rho*U) & mesh.Sf()
   //  );



    // Thermal conductivity field
    volScalarField lambda
    (
        IOobject
        (
            "lambda",
            runTime.timeName(),
            mesh
        ),
        alpha*lambdaL+(1.0-alpha)*lambdaS
    );

    // Heat capacity field
    volScalarField Cp
    (
        IOobject
        (
            "Cp",
            runTime.timeName(),
            mesh
        ),
       Foam::neg(T-Ts) * cpS + Foam::pos(T-Tl)*cpL+ Foam::pos(T-Ts)*Foam::neg(T-Tl)*hs/(Tl-Ts)
    );

    volScalarField diffusion_coeff
    (
        IOobject
        (
            "diffusion_coeff",
            runTime.timeName(),
            mesh
        ),
         lambda / Cp
    );

    // Kinematic viscosity field
    volScalarField nu
    (
        IOobject
        (
            "nu",
            runTime.timeName(),
            mesh
        ),
        alpha*nuL+(1.0-alpha)*nuS
    );

// kinematic turbulent thermal thermal conductivity m2/s
//Info<< "Reading field alphat\n" << endl;
//volScalarField alphat
//(
//    IOobject
 //   (
 //       "alphat",
 //       runTime.timeName(),
  //      mesh,
  //      IOobject::MUST_READ,
 //       IOobject::AUTO_WRITE
//    ),
//    mesh
//);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rhok*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p_rgh.name());

#include "createMRF.H"
#include "createIncompressibleRadiationModel.H"
#include "createFvOptions.H"
