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

#include "laminarFlux.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASSFModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(laminarFlux, 0);
addToRunTimeSelectionTable(RASSFModel, laminarFlux, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

laminarFlux::laminarFlux
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& C,
    const word& scalarFluxModelName,
    const word& modelName
)
:
    RASSFModel(modelName, U, phi, C, scalarFluxModelName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> laminarFlux::Dt() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Dt",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            //dimensionedScalar("Dt", nu()().dimensions(), 0.0)
            dimensionedScalar("Dt", dimensionSet(0,2,-1,0,0,0,0), 0.0)
        )
    );
}

tmp<volScalarField> laminarFlux::DEff() const
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

   // Info<<"Reading diffusivity "<< D << endl;

    //return tmp<volScalarField>(new volScalarField("DEff", D.value()));

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
            mesh_,
            dimensionedScalar("DEff", dimensionSet(0,2,-1,0,0,0,0), D.value())
        )
    );

}

/*tmp<fvVectorMatrix> laminar::divDevReff(volVectorField& U) const
{
    return
    (  
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}*/

tmp<fvScalarMatrix> laminarFlux::divDFlux(volScalarField& C) const
{
    return
    (
      -fvm::laplacian(DEff(), C) + fvm::SuSp(-fvc::div(phi()), C)
    );
}


bool laminarFlux::read()
{
    return RASSFModel::read();
}


void laminarFlux::correct()
{
    scalarFluxModel::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASSFModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
