/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Application
    interPhaseChangeFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids with phase-change
    (e.g. cavitation).  Uses a VOF (volume of fluid) phase-fraction based
    interface capturing approach, with optional mesh motion and mesh topology
    changes including adaptive re-meshing.

    The momentum and other fluid properties are of the "mixture" and a
    single momentum equation is solved.

    The set of phase-change models provided are designed to simulate cavitation
    but other mechanisms of phase-change are supported within this solver
    framework.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "subCycle.H"
//#include "./phaseChange/twoPhaseModels/interfaceProperties/interfaceProperties.H"
#include "./phaseChange/phaseChangeTwoPhaseMixtures/phaseChangeTwoPhaseMixture/phaseChangeTwoPhaseMixture.H"
#include "./phaseChange/smoothInterfaceProperties/smoothInterfaceProperties.H"
#include "kinematicMomentumTransportModel.H"
#include "pimpleControl.H"
//#include "pimpleMultiRegionControl.H" //multiRegionControl is only for dynamicFvMesh
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "regionProperties.H"

//solid
#include "solidThermo.H"
#include "./solid/solidRegionDiffNo.H"
#include "fluidThermo.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
//combined mesh generation
    #include "createMeshes.H"
    #include "createDyMControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

//Solid
    pimpleControl pimple2(mesh2);
//    pimpleMultiRegionControl pimples (mesh, mesh2);
    #include "./solid/createFields2.H"
//    #include "./solid/solidRegionDiffusionNo.H"	//아래 2개로 분배
    #include "./solid/readSolidTimeControls.H"

    turbulence->validate();

    #include "CourantNo.H"				//only for fluid
//For both meshes    
    #include "./solid/solidRegionDiffusionNo.H"		//1
    #include "./deltaT/setInitialMultiRegionDeltaT.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
      while (pimple2.run(runTime))
      {
 	#include "readDyMControls.H"
        #include "./solid/readSolidTimeControls.H"
        // Store divU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        volScalarField divU("divU0", fvc::div(fvc::absolute(phi, U)));

        #include "CourantNo.H"
        #include "alphaCourantNo.H"
	#include "./solid/solidRegionDiffusionNo.H"	//2
        #include "./deltaT/setMultiRegionDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    //mixture.correct();
		    interface.correct();
		    TPCmixture->correct();

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            divU = fvc::div(fvc::absolute(phi, U));

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            //mixture.correct();
	    //interface.correct();
	    TPCmixture->correct();

            #include "UEqn.H"			//Momentum-Predictor : get velocity vector
	    #include "EEqn.H"
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"		//Pressure-Corrector : get pressure gradient
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
	    //#include "EEqn.H"
	    //#include "./solid/setRegionSolidFields.H"
	    #include "./solid/solveSolid.H"     
	 
	}


        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
      }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
