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
    EPTT)LCR

Description
    Phan-Thin-Tanner (exponential form) non linear viscoelastic fluid model
    (Phan-Thien, N.; Tanner, R. I., 1977).

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    EPTT_LCR.C

\*---------------------------------------------------------------------------*/

#ifndef EPTT_LCR_H
#define EPTT_LCR_H

#include "viscoelasticLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscoelasticLaws
{

/*---------------------------------------------------------------------------*\
                           Class EPTT_LCR Declaration
\*---------------------------------------------------------------------------*/

class EPTT_LCR
:
    public viscoelasticLaw
{
    // Private data
    
        dictionary EPTT_LCR_Coeffs_;

        //- Transported viscoelastic stress
        volSymmTensorField tau_;

        //- LCR stress
        volSymmTensorField Psi;
    
        //- Eigen values
        volSymmTensorField E;

        //- Eigen vectors
        volTensorField R;

        // Model constants

        //- Density
        dimensionedScalar rho_;

        //- Solvent viscosity
        dimensionedScalar etaS_;

        //- Zero shear rate polymer viscosity
        dimensionedScalar etaP_;

        //- Extensibility parameter
        dimensionedScalar epsilon_;

        //- Relaxation time
        dimensionedScalar lambda_;

        //- Chains slip in relation to the continue environment
        dimensionedScalar zeta_;


public:

    //- Runtime type information
    TypeName("PTT-Exponential-LCR");

    // Constructors

        //- Construct from components
        EPTT_LCR
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

    virtual ~EPTT_LCR() = default;


    // Member Functions

        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const
        {
            return tau_;
        }

        //- Return the coupling term for the momentum equation
        virtual tmp<fvVectorMatrix> divTau(volVectorField& U) const;

        //- Correct the viscoelastic stress
        virtual void correct();
    
        //- Read viscoelasticProperties dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityLaws
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
