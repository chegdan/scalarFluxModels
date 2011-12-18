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

#include "scalarFluxModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(scalarFluxModel, 0);
defineRunTimeSelectionTable(scalarFluxModel, scalarFluxModel);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

scalarFluxModel::scalarFluxModel
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    //transportModel& transport,
    const volScalarField& C,
    const word& scalarFluxModelName
)
:
    regIOobject
    (
        IOobject
        (
            scalarFluxModelName,
            U.time().constant(),
            U.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    runTime_(U.time()),
    mesh_(U.mesh()),

    U_(U),
    phi_(phi),
   // transportModel_(transport),
    C_(C)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<scalarFluxModel> scalarFluxModel::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    //transportModel& transport,
    const volScalarField& C,
    const word& scalarFluxModelName
)
{
    // get model name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                "scalarTurbulenceProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("simulationType")
    );

    Info<< "Selecting scalar turbulence model type " << modelType << endl;

    scalarFluxModelConstructorTable::iterator cstrIter =
        scalarFluxModelConstructorTablePtr_->find(modelType);

    if (cstrIter == scalarFluxModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "scalarFluxModel::New(const volVectorField&, "
            "const surfaceScalarField&, transportModel&, const volScalarField&, const word&)"
        )   << "Unknown scalarFluxModel type "
            << modelType << nl << nl
            << "Valid scalarFluxModel types:" << endl
            << scalarFluxModelConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<scalarFluxModel>
    (
//        cstrIter()(U, phi, transport, C, scalarFluxModelName)
        cstrIter()(U, phi, C, scalarFluxModelName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void scalarFluxModel::correct()
{
//    transportModel_.correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
