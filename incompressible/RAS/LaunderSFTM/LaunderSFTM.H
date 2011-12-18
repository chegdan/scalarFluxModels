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
    Foam::incompressible::RASSFModels::LaunderSFTM

Description
    Uses the gradient diffusion hypothesis to model scalar flux in an 
    incompressible turbulent flow.

SourceFiles
    LaunderSFTM.C

\*---------------------------------------------------------------------------*/

#ifndef rasLaunderSFTM_H
#define rasLaunderSFTM_H

#include "RASSFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASSFModels
{

/*---------------------------------------------------------------------------*\
                           Class LaunderSFTM Declaration
\*---------------------------------------------------------------------------*/

class LaunderSFTM
:
    public RASSFModel
{

protected:

        // Model coefficients
            dimensionedScalar Sct_;

	    //dimensionedScalar D_

        // Fields
            volScalarField nut_;
            volScalarField Dt_;
//	    volScalarField k_;
//	    volScalarField epsilon_;
//	    volVectorField S_;
//	    volTensorField R_;

public:

    //- Runtime type information
    TypeName("LaunderSFTM");

    // Constructors

        //- Construct from components
        LaunderSFTM
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            //transportModel& transport,
	    const volScalarField& C,
            const word& scalarFluxModelName = scalarFluxModel::typeName,
            const word& modelName = typeName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<LaunderSFTM> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            //transportModel& transport,
	    const volScalarField& C,
            const word& scalarFluxModelName = scalarFluxModel::typeName
        );


    //- Destructor
    virtual ~LaunderSFTM()
    {}


    // Member Functions

        //- Return the turbulence diffusivity
        virtual tmp<volScalarField> Dt() const
	{
		return Dt_;
	}

        //- Return the effective diffusivity, i.e. the molecular diffusivity
        virtual tmp<volScalarField> DEff() const;

        //- Return the source term for the scalar equation
        virtual tmp<fvScalarMatrix> divDFlux(volScalarField& C) const;

        //- Correct the LaunderSFTM viscosity
        virtual void correct();

        //- Read turbulenceProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASSFModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
