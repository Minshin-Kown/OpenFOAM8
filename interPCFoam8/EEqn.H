Info<<"\n----------<Solve fluid hEqn>----------"<< endl;

{
	volScalarField& he = FluidThermo.he();
//	volScalarField he = C*T;



//---------------------------------22.01.13---------------------------------//

Pair<tmp<volScalarField> > vDotT = TPCmixture->vDotT();
//update volumetric condensation and vaporation rate to multiply (T-Tsat)									
		const volScalarField& vDotcT = vDotT[0]();
		const volScalarField& vDotvT = vDotT[1]();	


	volScalarField alphaDiffu = K/C;
//--------------------------------------------------------------------------//


{
    while (pimple.correctNonOrthogonal())
    {
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(rhoPhi, he) 		//phi is including rho [kg/s]
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
          -dpdt						//flow work term (Pdv devided by rho)
        )
      //+ thermophysicalTransport.divq(he)
        -fvm::laplacian (fvc::interpolate(alphaDiffu), he, "(alphaDiffu,he)")

     ==
       rho*(U&g)
      //+ rad.Sh(thermo, he)					//radiation term = off
      //+ reaction.Qdot()					//reaction term = off
      + fvOptions(rho, he)

    );

    EEqn.relax();

//--------------------------------22.01.13-----------------------------------//

    if (phaseChange){
		solve
		 (
			EEqn
			==
			//fvm::Sp(vDotvT-vDotcT,T) - (vDotvT-vDotcT)*TSat
			//fvm::Sp(vDotcT-vDotvT,T)	//----22.01.13 T가 없어서 불가능
			rho*fvm::Sp(vDotcT-vDotvT, he)
			- ((vDotcT-vDotvT)*TSatLocal)*C*rho
			
		 );
	 }
	 else
	 {
		 EEqn.solve();
	 }

//--------------------------------------------------------------------------//
    fvOptions.constrain(EEqn);

//    EEqn.solve();

    fvOptions.correct(he);
    }
}


    FluidThermo.correct();
    //rad.correct();

    Info<< "Min/max T:" << min(FluidThermo.T()).value() << ' '
        << max(FluidThermo.T()).value() << endl;
        
//------------------------------22.01.13---------------------------------------//

Pair<tmp<volScalarField> > mDot = TPCmixture->mDotT();
//update mass condensation and vaporation rate to multiply (T-Tsat)
	const volScalarField& mDotcT = mDot[0]();	//condensation (vapor)
	const volScalarField& mDotvT = mDot[1]();	//evaporation (fluid)
//		mDotT = (mDotcT-mDotvT)*(he-TSatLocal);
		mDotT = (mDotcT-mDotvT)*(FluidThermo.T()-TSatLocal);	//m dot fg = vapor - fluid

//-----------------------------------------------------------------------------//
Info<<"\n-----------------------------------------"<< endl;


