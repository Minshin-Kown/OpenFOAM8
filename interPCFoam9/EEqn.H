Info<<"\n----------<Solve fluid hEqn>----------"<< endl;

{
//	volScalarField& he = FluidThermo.he();




//---------------------------------22.01.13---------------------------------//

Pair<tmp<volScalarField> > vDotT = TPCmixture->vDotT();
//update volumetric condensation and vaporation rate to multiply (T-Tsat)									
		const volScalarField& vDotcT = vDotT[0]();
		const volScalarField& vDotvT = vDotT[1]();	


//	volScalarField alphaDiffu = K/C;
	volScalarField alphaDiffuT = K/(C*rho);
//--------------------------------------------------------------------------//


{
    while (pimple.correctNonOrthogonal())
    {
    fvScalarMatrix EEqn
    (
//        fvm::ddt(rho, he) + fvm::div(rhoPhi, he) 		//phi is including rho [kg/s]
	fvm::ddt(T) + fvm::div(phi, T)	


	//+ fvm::ddt(he) + fvm::div(phi, he)
      //+ fvc::ddt(rho, KE) + fvc::div(rhoPhi, KE)			//K is (1/2)u^2; K.E. term //he[J/kg]
      + (
         //   he.name() == "e"
         // ?
	 // fvc::div						//internal energy term
         //   (
         //       fvc::absolute(rhoPhi/fvc::interpolate(rho), U), 	
         //       p,
         //       "div(phiv,p)"
	 // )
          -dpdt/(rho*C)						//flow work term (Pdv devided by rho)
        )
//        -fvm::laplacian (fvc::interpolate(alphaDiffu), he, "(alphaDiffu,he)")
	-fvm::laplacian (fvc::interpolate(alphaDiffuT), T, "(alphaDiffu,he)")

     ==
       (U&g)/C
//      + fvOptions(rho, he)
        + fvOptions(T)

    );

    EEqn.relax();

//--------------------------------22.01.13-----------------------------------//

    if (phaseChange){
		 solve
		 (
			EEqn
			==
//			rho*fvm::Sp(vDotcT-vDotvT, he)
//			- ((vDotcT-vDotvT)*TSatLocal)*C*rho
			fvm::Sp(vDotcT-vDotvT, T)
			 - ((vDotcT-vDotvT)*TSatLocal)
			
		 );
	 }
	 else
	 {
		 EEqn.solve();
	 }

//--------------------------------------------------------------------------//
    fvOptions.constrain(EEqn);


    fvOptions.correct(T);
    }
}


//    FluidThermo.correct();	//------22.01.26 for T Eqn


    Info<< "Min/max T:" << min(T).value() << ' '
        << max(T).value() << endl;
        
//------------------------------22.01.13---------------------------------------//

Pair<tmp<volScalarField> > mDot = TPCmixture->mDotT();
//update mass condensation and vaporation rate to multiply (T-Tsat)
	const volScalarField& mDotcT = mDot[0]();	//condensation (vapor)
	const volScalarField& mDotvT = mDot[1]();	//evaporation (fluid)
//		mDotT = (mDotcT-mDotvT)*(he-TSatLocal);
		mDotT = (mDotcT-mDotvT)*((T)-TSatLocal);	//m dot fg = vapor - fluid

//-----------------------------------------------------------------------------//
Info<<"\n-----------------------------------------"<< endl;


