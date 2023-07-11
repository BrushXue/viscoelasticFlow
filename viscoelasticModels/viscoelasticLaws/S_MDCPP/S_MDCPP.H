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
    S-MDCPP

Description
    S-MDCPP non linear viscoelastic fluid model:
    Wei Wang, Xikui Li and Xianhong Han,
    Journal of Non-Newtonian Fluid Mechanics, 165,1480-1493, 2010.

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    S_MDCPP.C

\*---------------------------------------------------------------------------*/

#ifndef S_MDCP_H
#define S_MDCP_H

#include "viscoelasticLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscoelasticLaws
{

/*---------------------------------------------------------------------------*\
                           Class S_MDCPP Declaration
\*---------------------------------------------------------------------------*/

class S_MDCPP
:
    public viscoelasticLaw
{
    // Private data
        
        dictionary S_MDCPPCoeffs_;
        
        //- Transported viscoelastic stress
        volSymmTensorField tau_;

        // Model constants

        //- Density
        dimensionedScalar rho_;

        //- Solvent viscosity
        dimensionedScalar etaS_;

        //- Zero shear rate polymer viscosity
        dimensionedScalar etaP_;

        //- Material parameter, defining the amount of anisotropy
        dimensionedScalar zeta_;

        //- Relaxation time of the backbone tube orientation
        dimensionedScalar lambdaOb_;

        //- Relaxation time for the stretch
        dimensionedScalar lambdaOs_;

        //- Amount of arms at the end of a backbone
        dimensionedScalar q_;


public:

    //- Runtime type information
    TypeName("S_MDCPP");

    // Constructors

        //- Construct from components
        S_MDCPP
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~S_MDCPP() = default;


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
