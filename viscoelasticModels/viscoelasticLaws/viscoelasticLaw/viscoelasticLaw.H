/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    viscoelasticLaw

Description
    Abstract base class for viscoelastic flow models.

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    viscoelasticLaw.C
    newViscoelasticLaw.C

\*---------------------------------------------------------------------------*/

#ifndef viscoelasticLaw_H
#define viscoelasticLaw_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"

#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class viscoelasticLaw Declaration
\*---------------------------------------------------------------------------*/

class viscoelasticLaw
{
    
protected:
    
    // Protected data
    
    word name_;
    dictionary dict_;
    
    const volVectorField& U_;
    const surfaceScalarField& phi_;
    
    void decompose
    (
        const volTensorField& M,
        const volSymmTensorField& E, 
        const volTensorField& R,
        volSymmTensorField& S
    );

    void eigen
    (
        const volSymmTensorField& Psi_,
        volSymmTensorField& E,
        volTensorField& R
    );
    
    // Private Member Functions
    
    //- No copy construct
    viscoelasticLaw(const viscoelasticLaw&) = delete;
    
    //- No copy assignment
    void operator=(const viscoelasticLaw&) = delete;


public:

    //- Runtime type information
    TypeName("viscoelasticLaw");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            viscoelasticLaw,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (name, dict, U, phi)
        );

    // Selectors
    
    //- Return a reference to the selected viscoelastic law
    static autoPtr<viscoelasticLaw> New
    (
        const word& name,
        const dictionary& dict,
        const volVectorField& U,
        const surfaceScalarField& phi
    );
    

    // Constructors

        //- Construct from components
        viscoelasticLaw
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor
    virtual ~viscoelasticLaw() = default;
    

    // Member Functions

    //- Return the velocity field
    const volVectorField& U() const
    {
        return U_;
    }
    
    //- Return the face flux field
    const surfaceScalarField& phi() const
    {
        return phi_;
    }
    
        //- Return the phase transport properties dictionary
        const dictionary& dict() const
        {
            return dict_;
        }
    
        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const = 0;

        //- Return the coupling term for the momentum equation
        virtual tmp<fvVectorMatrix> divTau(volVectorField& U) const = 0;

        //- Correct the viscoelastic stress
        virtual void correct() = 0;
    
        //- Read viscoelasticProperties dictionary
        virtual bool read(const dictionary& dict) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
