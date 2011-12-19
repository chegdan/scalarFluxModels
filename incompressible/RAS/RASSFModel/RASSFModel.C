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

#include "RASSFModel.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(RASSFModel, 0);
defineRunTimeSelectionTable(RASSFModel, dictionary);
addToRunTimeSelectionTable(scalarFluxModel, RASSFModel, scalarFluxModel);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void RASSFModel::printCoeffs()
{
    if (printCoeffs_)
    {
        Info<< type() << "Coeffs" << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

RASSFModel::RASSFModel
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& C,
    const word& scalarFluxModelName
)
:
    scalarFluxModel(U, phi, C, scalarFluxModelName),

    IOdictionary
    (
        IOobject
        (
            "RASSFProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    turbulence_(lookup("turbulence")),
    activeScalar_(lookup("activeScalar")),
    printCoeffs_(lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(subOrEmptyDict(type + "Coeffs"))

   // kMin_("kMin", sqr(dimVelocity), SMALL),
    //epsilonMin_("epsilonMin", kMin_.dimensions()/dimTime, SMALL),
    //omegaMin_("omegaMin", dimless/dimTime, SMALL),

    //y_(mesh_)
{
    //kMin_.readIfPresent(*this);
    //epsilonMin_.readIfPresent(*this);
    //omegaMin_.readIfPresent(*this);

    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<RASSFModel> RASSFModel::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
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
                "RASSFProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("RASSFModel")
    );

    Info<< "Selecting RAS scalar flux model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "RASSFModel::New"
            "("
                "const volVectorField&, "
                "const surfaceScalarField&, "
                "const volScalarField&, "
                "const word&"
            ")"
        )   << "Unknown RASSFModel type "
            << modelType << nl << nl
            << "Valid RASSFModel types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<RASSFModel>
    (
        cstrIter()(U, phi, C, scalarFluxModelName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
scalar RASSFModel::yPlusLam(const scalar kappa, const scalar E) const
{
    scalar ypl = 11.0;

    for (int i=0; i<10; i++)
    {
        ypl = log(max(E*ypl, 1))/kappa;
    }

    return ypl;
}*/

void RASSFModel::correct()
{
    scalarFluxModel::correct();

 //   if (turbulence_ && mesh_.changing())
 //   {
 //       y_.correct();
 //   }
}


bool RASSFModel::read()
{
    //if (regIOobject::read())

    // Bit of trickery : we are both IOdictionary ('RASProperties') and
    // an regIOobject from the scalarFluxModel level. Problem is to distinguish
    // between the two - we only want to reread the IOdictionary.
    
    bool ok = IOdictionary::readData
    (
        IOdictionary::readStream
        (
            IOdictionary::type()
        )
    );
    IOdictionary::close();

    if (ok)
    {
        lookup("turbulence") >> turbulence_;

        if (const dictionary* dictPtr = subDictPtr(type() + "Coeffs"))
        {
            coeffDict_ <<= *dictPtr;
        }

      //  kMin_.readIfPresent(*this);
      //  epsilonMin_.readIfPresent(*this);
      //  omegaMin_.readIfPresent(*this);

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
