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
Author
    Dan Combest

Namespace
    Foam::incompressible::scalarFluxModels

Description
    Namespace for incompressible scalar flux models.

Class
    Foam::incompressible::scalarFluxModel

Description
    Abstract base class for incompressible scalar flux models
    (RAS, LES and laminar).

SourceFiles
    scalarFluxModel.C
    newscalarFluxModel.C

\*---------------------------------------------------------------------------*/

#ifndef scalarFluxModel_H
#define scalarFluxModel_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "incompressible/transportModel/transportModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

namespace incompressible
{

/*---------------------------------------------------------------------------*\
                           Class scalarFluxModel Declaration
\*---------------------------------------------------------------------------*/

class scalarFluxModel
:
    public regIOobject
{

protected:

    // Protected data

        const Time& runTime_;
        const fvMesh& mesh_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        //transportModel& transportModel_;

	//const dimensionedScalar& Dm_;

	const volScalarField& C_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        scalarFluxModel(const scalarFluxModel&);

        //- Disallow default bitwise assignment
        void operator=(const scalarFluxModel&);


public:

    //- Runtime type information
    TypeName("scalarFluxModel");


    // Declare run-time New selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            scalarFluxModel,
            scalarFluxModel,
            (
                const volVectorField& U,
                const surfaceScalarField& phi,
                //transportModel& transport,
		const volScalarField& C,
                const word& scalarFluxModelName
            ),
            //(U, phi, transport, C, scalarFluxModelName)
            (U, phi, C, scalarFluxModelName)
        );


    // Constructors

        //- Construct from components
        scalarFluxModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
           // transportModel& transport,
	    const volScalarField& C,
            const word& scalarFluxModelName = typeName

        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<scalarFluxModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            //transportModel& transport,
	    const volScalarField& C,
            const word& scalarFluxModelName = typeName
        );


    //- Destructor
    virtual ~scalarFluxModel()
    {}


    // Member Functions

        //- Access function to velocity field
        inline const volVectorField& U() const
        {
            return U_;
        }

        //- Access function to flux field
        inline const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- Access function to incompressible transport model
      /*  inline transportModel& transport() const
        {
            return transportModel_;
        }*/

        inline const volScalarField& C() const
        {
            return C_;
        }

        //- Return the laminar viscosity
        /*inline tmp<volScalarField> nu() const
        {
            return transportModel_.nu();
        }*/

	//- Return the turbulent diffusivity 
        virtual tmp<volScalarField> Dt() const = 0;

        //- Return the effective viscosity
        virtual tmp<volScalarField> DEff() const = 0;

        //- Return the source term for the scalar equation
        virtual tmp<fvScalarMatrix> divDFlux(volScalarField& C) const = 0;

/*
        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const = 0;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const = 0;

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const = 0;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const = 0;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const = 0;
*/

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct() = 0;

        //- Read LESProperties or RASProperties dictionary
        virtual bool read() = 0;

        //- Default dummy write function
        virtual bool writeData(Ostream&) const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
