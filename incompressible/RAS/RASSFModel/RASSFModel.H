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

Namespace
    Foam::incompressible::RASSFModels

Description
    Namespace for incompressible RAS turbulence models.

Class
    Foam::incompressible::RASSFModel

Description
    Abstract base class for incompressible turbulence models.

SourceFiles
    RASSFModel.C

\*---------------------------------------------------------------------------*/

#ifndef RASSFModel_H
#define RASSFModel_H

//#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "incompressible/scalarFluxModel/scalarFluxModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "nearWallDist.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
//#include "incompressible/transportModel/transportModel.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "bound.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
                           Class RASSFModel Declaration
\*---------------------------------------------------------------------------*/

class RASSFModel
:
    public scalarFluxModel,
    public IOdictionary
{

protected:

    // Protected data

        //- Turbulence on/off flag
        Switch turbulence_;

        //- Flag to print the model coeffs at run-time
        Switch printCoeffs_;

        //- Model coefficients dictionary
        dictionary coeffDict_;

        //- Lower limit of k
        //dimensionedScalar kMin_;

        //- Lower limit of epsilon
        //dimensionedScalar epsilonMin_;

        //- Lower limit for omega
        //dimensionedScalar omegaMin_;

        //- Near wall distance boundary field
        //nearWallDist y_;


    // Protected Member Functions

        //- Print model coefficients
        virtual void printCoeffs();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        RASSFModel(const RASSFModel&);

        //- Disallow default bitwise assignment
        void operator=(const RASSFModel&);


public:

    //- Runtime type information
    TypeName("RASSFModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            RASSFModel,
            dictionary,
            (
                const volVectorField& U,
                const surfaceScalarField& phi,
		const volScalarField& C,
                const word& scalarFluxModelName
            ),
            (U, phi, C, scalarFluxModelName)
        );


    // Constructors

        //- Construct from components
        RASSFModel
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi,
	    const volScalarField& C,
            const word& scalarFluxModelName = scalarFluxModel::typeName

        );


    // Selectors

        //- Return a reference to the selected RAS model
        static autoPtr<RASSFModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
	    const volScalarField& C,
            const word& scalarFluxModelName = scalarFluxModel::typeName

        );


    //- Destructor
    virtual ~RASSFModel()
    {}


    // Member Functions

        // Access

     /*       //- Return the lower allowable limit for k (default: SMALL)
            const dimensionedScalar& kMin() const
            {
                return kMin_;
            }

            //- Return the lower allowable limit for epsilon (default: SMALL)
            const dimensionedScalar& epsilonMin() const
            {
                return epsilonMin_;
            }

            //- Return the lower allowable limit for omega (default: SMALL)
            const dimensionedScalar& omegaMin() const
            {
                return omegaMin_;
            }

            //- Allow kMin to be changed
            dimensionedScalar& kMin()
            {
                return kMin_;
            }

            //- Allow epsilonMin to be changed
            dimensionedScalar& epsilonMin()
            {
                return epsilonMin_;
            }

            //- Allow omegaMin to be changed
            dimensionedScalar& omegaMin()
            {
                return omegaMin_;
            }

            //- Return the near wall distances
            const nearWallDist& y() const
            {
                return y_;
            }

            //- Calculate y+ at the edge of the laminar sublayer
            //scalar yPlusLam(const scalar kappa, const scalar E) const;
*/
            //- Const access to the coefficients dictionary
            const dictionary& coeffDict() const
            {
                return coeffDict_;
            }


        //- Return the effective viscosity
        virtual tmp<volScalarField> DEff() const
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
                new volScalarField("DEff", Dt() + D.value())
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
