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

Class
    Foam::incompressible::laminarFlux

Description
    Scalar Flux model for laminarFlux incompressible flow.

SourceFiles
    laminarFlux.C

\*---------------------------------------------------------------------------*/

#ifndef laminarFlux_H
#define laminarFlux_H

#include "scalarFluxModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
                           Class laminarFlux Declaration
\*---------------------------------------------------------------------------*/

class laminarFlux
:
    public scalarFluxModel
{

public:

    //- Runtime type information
    TypeName("laminarFlux");

    // Constructors

        //- Construct from components
        laminarFlux
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            //transportModel& transport,
	    const volScalarField& C,
            const word& scalarFluxModelName = scalarFluxModel::typeName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<laminarFlux> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            //transportModel& transport,
	    const volScalarField& C,
            const word& scalarFluxModelName = scalarFluxModel::typeName
        );


    //- Destructor
    virtual ~laminarFlux()
    {}


    // Member Functions

        //- Return the turbulence diffusivity, i.e. 0 for laminarFlux flow
        virtual tmp<volScalarField> Dt() const;

        //- Return the effective diffusivity, i.e. the molecular diffusivity
        virtual tmp<volScalarField> DEff() const;

        //- Return the source term for the scalar equation
        virtual tmp<fvScalarMatrix> divDFlux(volScalarField& C) const;

        //- Correct the laminarFlux viscosity
        virtual void correct();

        //- Read turbulenceProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
