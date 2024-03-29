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
    XPP_SE

Description
    eXtended Pom-Pom (XPP) single equation non linear viscoelastic fluid model
    (Verbeeten, W. M. H., Computational Polymer Melt Rheology, 2001).

Author
    Jovani L. Favero.  All rights reserved

SourceFiles
    XPP_SE.C

\*---------------------------------------------------------------------------*/

#ifndef XPP_SE_H
#define XPP_SE_H

#include "viscoelasticLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscoelasticLaws
{

/*---------------------------------------------------------------------------*\
                           Class XPP_SE Declaration
\*---------------------------------------------------------------------------*/

class XPP_SE
:
    public viscoelasticLaw
{
    // Private data

        dictionary XPP_SECoeffs_;
    
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
        dimensionedScalar alpha_;

        //- Relaxation time of the backbone tube orientation
        dimensionedScalar lambdaOb_;

        //- Relaxation time for the stretch
        dimensionedScalar lambdaOs_;

        //- Amount of arms at the end of a backbone
        dimensionedScalar q_;


public:

    //- Runtime type information
    TypeName("XPP-SE");

    // Constructors

        //- Construct from components
        XPP_SE
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~XPP_SE() = default;


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
