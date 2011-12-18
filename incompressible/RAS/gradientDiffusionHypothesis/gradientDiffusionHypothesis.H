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
    Foam::incompressible::RASSFModels::gradientDiffusionHypothesis

Description
    Uses the gradient diffusion hypothesis to model scalar flux in an 
    incompressible turbulent flow.

SourceFiles
    gradientDiffusionHypothesis.C

\*---------------------------------------------------------------------------*/

#ifndef rasgradientDiffusionHypothesis_H
#define rasgradientDiffusionHypothesis_H

#include "RASSFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASSFModels
{

/*---------------------------------------------------------------------------*\
                           Class gradientDiffusionHypothesis Declaration
\*---------------------------------------------------------------------------*/

class gradientDiffusionHypothesis
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

public:

    //- Runtime type information
    TypeName("gradientDiffusionHypothesis");

    // Constructors

        //- Construct from components
        gradientDiffusionHypothesis
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
        static autoPtr<gradientDiffusionHypothesis> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            //transportModel& transport,
	    const volScalarField& C,
            const word& scalarFluxModelName = scalarFluxModel::typeName
        );


    //- Destructor
    virtual ~gradientDiffusionHypothesis()
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

        //- Correct the gradientDiffusionHypothesis viscosity
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
