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

\*---------------------------------------------------------------------------*/

#include "Oldroyd_B.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscoelasticLaws
    {
        defineTypeNameAndDebug(Oldroyd_B, 0);
        addToRunTimeSelectionTable(viscoelasticLaw, Oldroyd_B, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscoelasticLaws::Oldroyd_B::Oldroyd_B
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscoelasticLaw(name, dict, U, phi),

    Oldroyd_BCoeffs_
    (
        dict.optionalSubDict(typeName + "Coeffs")
    ),

    tau_
    (
        IOobject
        (
            "tau" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    rho_("rho", dimDensity, Oldroyd_BCoeffs_),
    etaS_("etaS", dimDynamicViscosity, Oldroyd_BCoeffs_),
    etaP_("etaP", dimDynamicViscosity, Oldroyd_BCoeffs_),
    lambda_("lambda", dimTime, Oldroyd_BCoeffs_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::viscoelasticLaws::Oldroyd_B::divTau(volVectorField& U) const
{
    return
    (
        fvc::div(tau_ / rho_, "div(tau)")
      - fvc::div(etaP_ / rho_ * fvc::grad(U), "div(grad(U))")
      + fvm::laplacian( (etaP_ + etaS_)/rho_, U, "laplacian(eta,U)")
    );
}


void Foam::viscoelasticLaws::Oldroyd_B::correct()
{
    // Velocity gradient tensor
    volTensorField L(fvc::grad(U()));
    
    // Convected derivate term
    volTensorField C(tau_ & L);

    // Twice the rate of deformation tensor
    volSymmTensorField twoD(twoSymm(L));

    // Stress transport equation
    fvSymmTensorMatrix tauEqn
    (
        fvm::ddt(tau_)
      + fvm::div(phi(), tau_)
     ==
        etaP_/lambda_*twoD
      + twoSymm(C)
      - fvm::Sp(1/lambda_, tau_)
    );

    tauEqn.relax();
    tauEqn.solve();
}

bool Foam::viscoelasticLaws::Oldroyd_B::read
(
    const dictionary& dict
)
{
    viscoelasticLaw::read(dict);
    
    Oldroyd_BCoeffs_ = dict.optionalSubDict(typeName + "Coeffs");
    
    Oldroyd_BCoeffs_.readEntry("rho", rho_);
    Oldroyd_BCoeffs_.readEntry("etaS", etaS_);
    Oldroyd_BCoeffs_.readEntry("etaP", etaP_);
    Oldroyd_BCoeffs_.readEntry("lambda", lambda_);
    
    return true;
}
// ************************************************************************* //
