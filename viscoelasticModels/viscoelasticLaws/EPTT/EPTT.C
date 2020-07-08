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

#include "EPTT.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscoelasticLaws
    {
        defineTypeNameAndDebug(EPTT, 0);
        addToRunTimeSelectionTable(viscoelasticLaw, EPTT, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscoelasticLaws::EPTT::EPTT
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscoelasticLaw(name, dict, U, phi),

    EPTTCoeffs_
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

    rho_("rho", dimDensity, EPTTCoeffs_),
    etaS_("etaS", dimDynamicViscosity, EPTTCoeffs_),
    etaP_("etaP", dimDynamicViscosity, EPTTCoeffs_),
    epsilon_("epsilon", dimless, EPTTCoeffs_),
    lambda_("lambda", dimTime, EPTTCoeffs_),
    zeta_("zeta", dimless, EPTTCoeffs_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::viscoelasticLaws::EPTT::divTau(volVectorField& U) const
{
    dimensionedScalar etaPEff = etaP_;

    return
    (
        fvc::div(tau_/rho_, "div(tau)")
      - fvc::laplacian(etaPEff/rho_, U, "laplacian(etaPEff,U)")
      + fvm::laplacian( (etaPEff + etaS_)/rho_, U, "laplacian(etaPEff+etaS,U)")
    );
}


void Foam::viscoelasticLaws::EPTT::correct()
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
      - zeta_*symm(tau_ & twoD)
      - fvm::Sp
        (
            (1/lambda_)*Foam::exp(epsilon_*lambda_/etaP_*tr(tau_)),
            tau_
        )
    );

    tauEqn.relax();
    tauEqn.solve();
}

bool Foam::viscoelasticLaws::EPTT::read
(
    const dictionary& dict
)
{
    viscoelasticLaw::read(dict);
    
    EPTTCoeffs_ = dict.optionalSubDict(typeName + "Coeffs");
    
    EPTTCoeffs_.readEntry("rho", rho_);
    EPTTCoeffs_.readEntry("etaS", etaS_);
    EPTTCoeffs_.readEntry("etaP", etaP_);
    EPTTCoeffs_.readEntry("epsilon", epsilon_);
    EPTTCoeffs_.readEntry("lambda", lambda_);
    EPTTCoeffs_.readEntry("zeta", zeta_);
    
    return true;
}
// ************************************************************************* //
