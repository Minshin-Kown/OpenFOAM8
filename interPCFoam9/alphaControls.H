//-----------------------22.01.25---------------------------------------------------//
  
  Pair<tmp<volScalarField> > mDot = TPCmixture->mDotAlphal();
				const volScalarField& mDotcAlpha= mDot[0]();    //Vapor
				const volScalarField& mDotvAlpha= mDot[1]();    //Liquid

				mDotAlpha = mDotvAlpha*alpha1+mDotcAlpha*(1.0-alpha1);
                                    //alpha1=1이면 liquid
//-----------------------------------------------------------------------------------//
/*
surfaceScalarField rhoPhi
     (
        IOobject
           (
             "rhoPhi",
             runTime.timeName(),
             mesh
           ),
         mesh,
         dimensionedScalar(dimMass/dimTime, 0)
     );
*/

const dictionary& alphaControls = mesh.solverDict(alpha1.name());	//fvSolution에서 alpha1.name 읽어와 별칭 alphaControls

const label nAlphaCorr(alphaControls.lookup<label>("nAlphaCorr"));	//alphaControls의 nAlphaCoor

const label nAlphaSubCycles(alphaControls.lookup<label>("nAlphaSubCycles"));

//const bool MULESCorr(alphaControls.lookupOrDefault<Switch>("MULESCorr", false));	//----trun off 22.01.14

// Apply the compression correction from the previous iteration
// Improves efficiency for steady-simulations but can only be applied
// once the alpha field is reasonably steady, i.e. fully developed
// const bool alphaApplyPrevCorr
// (
//     alphaControls.lookupOrDefault<Switch>("alphaApplyPrevCorr", false)
// );
