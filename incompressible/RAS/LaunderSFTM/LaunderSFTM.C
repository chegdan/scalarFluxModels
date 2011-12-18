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

\*---------------------------------------------------------------------------*/

#include "LaunderSFTM.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASSFModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LaunderSFTM, 0);
addToRunTimeSelectionTable(RASSFModel, LaunderSFTM, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LaunderSFTM::LaunderSFTM
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& C,
    const word& scalarFluxModelName,
    const word& modelName
)
:
    RASSFModel(modelName, U, phi, C, scalarFluxModelName),

    Sct_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sct",
            coeffDict_,
            1
        )
    ),

    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    Dt_
    (
        IOobject
        (
            "Dt",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
         mesh_
    )

{
    read();

    Dt_ = nut_/Sct_;//turbulent diffusivity

    Dt_.correctBoundaryConditions();
    
    printCoeffs();

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> LaunderSFTM::DEff() const
{
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime_.constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar D
    (
        transportProperties.lookup("D")
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "DEff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
	    D + Dt()
        )
    );

}

tmp<fvScalarMatrix> LaunderSFTM::divDFlux(volScalarField& C) const
{
    return
    (
      -fvm::laplacian(DEff(), C) + fvm::SuSp(-fvc::div(phi()), C)
    );
}


bool LaunderSFTM::read()
{
    if (RASSFModel::read())
    {
	//Info<<"Reading turbulent Schmidt number"<< endl;

        Sct_.readIfPresent(coeffDict());

	//Info<<"Using  turbulent Schmidt of "<< Sct_<< endl;

        return true;
    }
    else
    {
        return false;
    }
}


void LaunderSFTM::correct()
{
    RASSFModel::correct();

    if (!turbulence_)
    {
	Info<<"turbulence turned off"<< endl;
        return;
    }

    //recalculate the turbulent diffusivity
    Dt_ = nut_/Sct_;
    Dt_.correctBoundaryConditions();

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASSFModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
