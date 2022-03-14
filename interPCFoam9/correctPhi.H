/*
CorrectPhi
(
    U,
    phi,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU())),
    divU,
    pimple,
    true
);

#include "continuityErrs.H"
*/

{
    #include "continuityErrs.H"

    wordList pcorrTypes                             //인스턴스 선언 및 정의 int b(10); 꼴 @C++
    (
        p_rgh.boundaryField().size(),           
        zeroGradientFvPatchScalarField::typeName
    );

    forAll (p_rgh.boundaryField(), i)                //앞서 선언된 i개의 boundary field에 대해 모두
    {
        if (p_rgh.boundaryField()[i].fixesValue())      //fixesValue라면, 하기와 같은 typeName을 저장
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr                               //volScalarField의 pcorr 인스턴스 (장) 정의
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", p_rgh.dimensions(), 0.0),
        pcorrTypes
    );

    dimensionedScalar rAUf("(1|A(U))", dimTime/rho.dimensions(), 1.0);  //rAUf 인스턴스 선언 및 정의

    adjustPhi(phi, U, pcorr);   

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn                     //pcorrEqn:pressure correction due to mesh movement. 
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi)
        );

        pcorrEqn.setReference(pRefCell, pRefValue);
        pcorrEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }
    }

    #include "continuityErrs.H"
}
