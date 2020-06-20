{



    fvScalarMatrix hEqn
    
    (
	fvm::ddt(rho,enthalpy)
	+fvm::div(phi*fvc::interpolate(rho),enthalpy)
	-fvm::laplacian(diffusion_coeff,enthalpy)
    );


    hEqn.relax();
    hEqn.solve();

  //  enthalpy.correctBoundaryConditions();

      fvOptions.correct(enthalpy) ;

      dimensionedScalar      criteria_1(   cpS*Ts );

      dimensionedScalar      criteria_2(   cpS*Ts + hs  );
      
   // Calculate the temperature field
      T= Foam::pos(criteria_1-enthalpy)*(enthalpy/cpS)
      	+Foam::pos(enthalpy-criteria_2)*(Tl+(enthalpy-criteria_2)/cpL)
      	+Foam::pos(enthalpy-criteria_1)*Foam::pos(criteria_2-enthalpy)*(Ts+ (enthalpy-criteria_1)*(Tl-Ts)/hs);

      T.correctBoundaryConditions();

    // Calculate the phase field
      alpha = max(    min((T-Ts)/(Tl-Ts),scalar(1))   ,scalar(0));


      scalar tv = runTime.value();

      alphaValues.insert(tv, alpha.weightedAverage(mesh.V()).value());

      Info<< "Total liquid fraction [-] = " << alpha.weightedAverage(mesh.V()).value()<< endl;


   // alpha=(T-Ts)/(Tl-Ts);
   // alpha= max(scalar(0.0),alpha);
  //  alpha= min(scalar(1.0),alpha);

    rho = alpha*rhoL+(1.0-alpha)*rhoS;

  // Calculate the mass conservative fluxes
  //surfaceScalarField rhophi("rhophi", linearInterpolate(rho*U) & mesh.Sf());

    lambda = alpha*lambdaL+(1.0-alpha)*lambdaS;

    Cp=Foam::neg(T-Ts)*cpS + Foam::pos(T-Tl)*cpL+ Foam::pos(T-Ts)*Foam::neg(T-Tl)*hs/(Tl-Ts);

    Cp.correctBoundaryConditions();


  //  diffusion_coeff = lambda / (Cp*rho);

    diffusion_coeff = lambda / Cp;

    nu = alpha*nuL+(1.0-alpha)*nuS;

 //   rhok = 1.0 - Foam::max(beta*(T - Tl),0.0);

    rhok = 1.0 - Foam::max(beta*(T - Ts),0.0);

    DC = DCl*Foam::pow(1.0-alpha,2)/(Foam::pow(alpha,3)+DCs);

    Info<< "min(T) " << gMin(T) << endl;

    Info<< "max(T) " << gMax(T) << endl;

    Info<< "min(enthalpy) " << gMin(enthalpy) << endl;

    Info<< "max(enthalpy) " << gMax(enthalpy) << endl;

}


