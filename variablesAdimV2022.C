/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    ProducionUU
Description
    Calculates and writes the uu velocity fluctiations production term
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvc.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

    // new
    #include "createControl.H"
    #include "createTimeControls.H"
    

	forAll(timeDirs, timeI)
      {
      runTime.setTime(timeDirs[timeI], timeI);

      Info<< "Time = " << runTime.timeName() << endl;

      IOobject dissipationTKEMeanHeader
      (
            "dissipationTKEMean",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
      );

      IOobject prodTKEMeanHeader
      (
            "prodTKEMean",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
      );
      
      IOobject Lambda2Header
      (
            "Lambda2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
      );      
      
      IOobject pressureDiffusionTKEMeanHeader
      (
            "pressureDiffusionTKEMean",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
      );      
      
      IOobject pressureStrainTKEMeanHeader
      (
            "pressureStrainTKEMean",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
      );        
           
      IOobject turbulenceTransportTKEMeanHeader
      (
            "turbulenceTransportTKEMean",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
      );      
      
      IOobject UHeader
      (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
      );      
      
      IOobject UMeanHeader
      (
            "UMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
      );      
      
      IOobject viscousDiffusionTKEMeanHeader
      (
            "viscousDiffusionTKEMean",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
      );      
      
      IOobject vorticityMeanHeader
      (
            "vorticityMean",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
      );      
      
      IOobject RMeanHeader
      (
            "RMean",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
      );

      IOobject UPrime2MeanHeader
      (
            "UPrime2MeanHeader",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
      );

      IOobject rhoHeader
      (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
      );
      
      IOobject muLamHeader
      (
            "muLam",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
      );
      
      IOobject muTHeader
      (
            "muT",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
      );
            
      IOdictionary transportProperties
      (
            IOobject
            (
                  "transportProperties",
                  runTime.constant(),
                  mesh,
                  IOobject::MUST_READ_IF_MODIFIED,
                  IOobject::NO_WRITE,
                  false
            )
      );

      dimensionedScalar nu
      (
          "nu",
          dimViscosity,
          transportProperties.lookup("nu")
      );


      if (dissipationTKEMeanHeader.typeHeaderOk<volScalarField>(true))
      {
            mesh.readUpdate();

            Info<< "    Reading dissipationTKEMean" << endl;
            volScalarField dissipationTKEMean(dissipationTKEMeanHeader, mesh);
            Info<< "    Calculating dissipationTKEMeanAdimUTau" << endl;

            volScalarField dissipationTKEMeanAdimUTau
            (
                  IOobject
                  (
                  "dissipationTKEMeanAdimUTau",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ
                  ),
				  
				  dissipationTKEMean/(pow(0.654476306,4)/nu)
            );
            dissipationTKEMeanAdimUTau.write();
      }
      else
      {
            Info<< "    No dissipationTKEMean" << endl;
      }


      if (prodTKEMeanHeader.typeHeaderOk<volScalarField>(true))
      {
            mesh.readUpdate();

            Info<< "    Reading prodTKEMean" << endl;
            volScalarField prodTKEMean(prodTKEMeanHeader, mesh);
            Info<< "    Calculating prodTKEMeanAdimUTau" << endl;

            volScalarField prodTKEMeanAdimUTau
            (
                  IOobject
                  (
                  "prodTKEMeanAdimUTau",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ
                  ),
                  prodTKEMean/(pow(0.654476306,4)/nu)
            );
            prodTKEMeanAdimUTau.write();
      }
      else
      {
            Info<< "    No dissipationTKEMean" << endl;
      }

      if (pressureDiffusionTKEMeanHeader.typeHeaderOk<volScalarField>(true))
      {
            mesh.readUpdate();

            Info<< "    Reading pressureDiffusionTKEMean" << endl;
            volScalarField pressureDiffusionTKEMean(pressureDiffusionTKEMeanHeader, mesh);
            Info<< "    Calculating pressureDiffusionTKEMeanAdimUTau" << endl;

            volScalarField pressureDiffusionTKEMeanAdimUTau
            (
                  IOobject
                  (
                  "pressureDiffusionTKEMeanAdimUTau",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ
                  ),
                  pressureDiffusionTKEMean/(pow(0.654476306,4)/nu)
            );
            pressureDiffusionTKEMeanAdimUTau.write();
      }
      else
      {
            Info<< "    No pressureDiffusionTKEMean" << endl;
      }

      if (pressureStrainTKEMeanHeader.typeHeaderOk<volScalarField>(true))
      {
            mesh.readUpdate();

            Info<< "    Reading pressureStrainTKEMean" << endl;
            volScalarField pressureStrainTKEMean(pressureStrainTKEMeanHeader, mesh);
            Info<< "    Calculating pressureStrainTKEMeanAdimUTau" << endl;

            volScalarField pressureStrainTKEMeanAdimUTau
            (
                  IOobject
                  (
                  "pressureStrainTKEMeanAdimUTau",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ
                  ),
                  pressureStrainTKEMean/(pow(0.654476306,4)/nu)
            );
            pressureStrainTKEMeanAdimUTau.write();
      }
      else
      {
            Info<< "    No pressureStrainTKEMean" << endl;
      }

      if (prodTKEMeanHeader.typeHeaderOk<volScalarField>(true))
      {
            mesh.readUpdate();

            Info<< "    Reading prodTKEMean" << endl;
            volScalarField prodTKEMean(prodTKEMeanHeader, mesh);
            Info<< "    Calculating prodTKEMeanAdimUTau" << endl;

            volScalarField prodTKEMeanAdimUTau
            (
                  IOobject
                  (
                  "prodTKEMeanAdimUTau",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ
                  ),
                  prodTKEMean/(pow(0.654476306,4)/nu)
            );
            prodTKEMeanAdimUTau.write();
      }
      else
      {
            Info<< "    No prodTKEMean" << endl;
      }

      if (turbulenceTransportTKEMeanHeader.typeHeaderOk<volScalarField>(true))
      {
            mesh.readUpdate();

            Info<< "    Reading turbulenceTransportTKEMean" << endl;
            volScalarField turbulenceTransportTKEMean(turbulenceTransportTKEMeanHeader, mesh);
            Info<< "    Calculating turbulenceTransportTKEMeanAdimUTau" << endl;

            volScalarField turbulenceTransportTKEMeanAdimUTau
            (
                  IOobject
                  (
                  "turbulenceTransportTKEMeanAdimUTau",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ
                  ),
                  turbulenceTransportTKEMean/(pow(0.654476306,4)/nu)
            );
            turbulenceTransportTKEMeanAdimUTau.write();
      }
      else
      {
            Info<< "    No turbulenceTransportTKEMean" << endl;
      }

      if (viscousDiffusionTKEMeanHeader.typeHeaderOk<volScalarField>(true))
      {
            mesh.readUpdate();

            Info<< "    Reading viscousDiffusionTKEMean" << endl;
            volScalarField viscousDiffusionTKEMean(viscousDiffusionTKEMeanHeader, mesh);
            Info<< "    Calculating viscousDiffusionTKEMeanAdimUTau" << endl;

            volScalarField viscousDiffusionTKEMeanAdimUTau
            (
                  IOobject
                  (
                  "viscousDiffusionTKEMeanAdimUTau",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ
                  ),
                  viscousDiffusionTKEMean/(pow(0.654476306,4)/nu)
            );
            viscousDiffusionTKEMeanAdimUTau.write();
      }
      else
      {
            Info<< "    No viscousDiffusionTKEMean" << endl;
      }

      if (UMeanHeader.typeHeaderOk<volVectorField>(true))
      {
            mesh.readUpdate();

            Info<< "    Reading UMean" << endl;
            volVectorField UMean(UMeanHeader, mesh);
                        
            Info<< "    Calculating UAdimUBulk" << endl;

            volVectorField UMeanAdimUBulk
            (
                  IOobject
                  (
                  "UMeanAdimUBulk",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ
                  ),
                  UMean/9.6144
            );
            UMeanAdimUBulk.write();
      }
      else
      {
            Info<< "    No UMean" << endl;
      }

      if (RMeanHeader.typeHeaderOk<volSymmTensorField>(true))
      {
            mesh.readUpdate();

            Info<< "    Reading RMean" << endl;
            volSymmTensorField RMean(RMeanHeader, mesh);
                        
            Info<< "    Calculating kAdimUBulk" << endl;

            volScalarField kAdimUBulk
            (
                  IOobject
                  (
                  "kAdimUBulk",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ
                  ),
                  0.5*tr(RMean)/pow(9.6144,2)
            );
            kAdimUBulk.write();
            
            Info<< "    Calculating turbulenceIntensity" << endl;
            volScalarField turbulenceIntensity
            (
                  IOobject
                  (
                  "turbulenceIntensity",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ
                  ),
                  pow((2/3*0.5*tr(RMean)),0.5)/9.6144
            );
            turbulenceIntensity.write();
			
            Info<< "    Calculating URMSAdimUBulk " << endl;
			volScalarField RMeanxx(RMean.component(symmTensor::XX));
			volScalarField RMeanyy(RMean.component(symmTensor::YY));			
			volScalarField RMeanzz(RMean.component(symmTensor::ZZ));
			
            volVectorField URMSAdimUBulk
            (
				IOobject
				(
				"URMSAdimUBulk",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ
				),
				mesh,			
				dimensionedVector
				(
					"URMSAdimUBulk", 
					dimVelocity, 
					vector(0.0,0.0,0.0)
				)	  
            );
			URMSAdimUBulk.component(0) =  pow(RMeanxx,0.5)/9.6144;
			URMSAdimUBulk.component(1) =  pow(RMeanyy,0.5)/9.6144;
			URMSAdimUBulk.component(2) =  pow(RMeanzz,0.5)/9.6144;			
			
            URMSAdimUBulk.write();						
      }
      else
      {
            Info<< "    No RMean" << endl;
      }
	  
      if (UPrime2MeanHeader.typeHeaderOk<volSymmTensorField>(true))
      {
            mesh.readUpdate();

            Info<< "    Reading UPrime2Mean" << endl;
            volSymmTensorField UPrime2Mean(UPrime2MeanHeader, mesh);
                        
            Info<< "    Calculating kAdimUBulk" << endl;

            volScalarField kAdimUBulk
            (
                  IOobject
                  (
                  "kAdimUBulk",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ
                  ),
                  0.5*tr(UPrime2Mean)/pow(9.6144,2)
            );
            kAdimUBulk.write();
            
            Info<< "    Calculating turbulenceIntensity" << endl;
            volScalarField turbulenceIntensity
            (
                  IOobject
                  (
                  "turbulenceIntensity",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ
                  ),
                  pow((2/3*0.5*tr(UPrime2Mean)),0.5)/9.6144
            );
            turbulenceIntensity.write();
			
            Info<< "    Calculating URMSAdimUBulk " << endl;
			volScalarField UPrime2Meanxx(UPrime2Mean.component(symmTensor::XX));
			volScalarField UPrime2Meanyy(UPrime2Mean.component(symmTensor::YY));			
			volScalarField UPrime2Meanzz(UPrime2Mean.component(symmTensor::ZZ));
			
            volVectorField URMSAdimUBulk
            (
				IOobject
				(
				"URMSAdimUBulk",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ
				),
				mesh,			
				dimensionedVector
				(
					"URMSAdimUBulk", 
					dimVelocity, 
					vector(0.0,0.0,0.0)
				)	  
            );
			URMSAdimUBulk.component(0) =  pow(UPrime2Meanxx,0.5)/9.6144;
			URMSAdimUBulk.component(1) =  pow(UPrime2Meanyy,0.5)/9.6144;
			URMSAdimUBulk.component(2) =  pow(UPrime2Meanzz,0.5)/9.6144;			
			
            URMSAdimUBulk.write();						
      }
      else
      {
            Info<< "    No UPrime2Mean" << endl;
      }
	  
      Info<< endl;
    }

return 0;
}


// ************************************************************************* //
