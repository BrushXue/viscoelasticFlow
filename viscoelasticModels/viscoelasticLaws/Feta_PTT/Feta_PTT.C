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

#include "Feta_PTT.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscoelasticLaws
    {
        defineTypeNameAndDebug(Feta_PTT, 0);
        addToRunTimeSelectionTable(viscoelasticLaw, Feta_PTT, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscoelasticLaws::Feta_PTT::Feta_PTT
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscoelasticLaw(name, dict, U, phi),

    Feta_PTTCoeffs_
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
    
    rho_("rho", dimDensity, Feta_PTTCoeffs_),
    etaS_("etaS", dimDynamicViscosity, Feta_PTTCoeffs_),
    etaP_("etaP", dimDynamicViscosity, Feta_PTTCoeffs_),
    epsilon_("epsilon", dimless, Feta_PTTCoeffs_),
    lambda_("lambda", dimTime, Feta_PTTCoeffs_),
    zeta_("zeta", dimless, Feta_PTTCoeffs_),
    A_("A", dimless, Feta_PTTCoeffs_),
    a_("a", dimless, Feta_PTTCoeffs_),
    b_("b", dimless, Feta_PTTCoeffs_),

    etaPEff_
    (
        IOobject
        (
            "etaPEff",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        etaP_/
        (
            Foam::pow(scalar(1.0) + A_*Foam::pow(0.5*(Foam::sqr(tr(tau_))
          - tr(tau_ & tau_))*Foam::sqr(lambda_)/Foam::sqr(etaP_), a_), b_)
        )
    ),

    lambdaEff_
    (
        IOobject
        (
            "lambdaEff",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        lambda_ / (scalar(1)  + epsilon_*lambda_*tr(tau_) / etaP_)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::viscoelasticLaws::Feta_PTT::divTau(volVectorField& U) const
{
    return
    (
        fvc::div(tau_ / rho_, "div(tau)")
      - fvc::div(etaP_ / rho_ * fvc::grad(U), "div(grad(U))")
      + fvm::laplacian((etaP_ + etaS_) / rho_, U, "laplacian(eta,U)")
    );
}


void Foam::viscoelasticLaws::Feta_PTT::correct()
{
    // Velocity gradient tensor
    volTensorField L(fvc::grad(U()));

    // Convected derivate term
    volTensorField C(tau_ & L);

    // Twice the rate of deformation tensor
    volSymmTensorField twoD(twoSymm(L));

    // etaP effective
    etaPEff_ = etaP_/
        (
            Foam::pow(scalar(1.0) + A_*Foam::pow(0.5*( Foam::sqr(tr(tau_))
          - tr(tau_ & tau_)) * Foam::sqr(lambda_)/Foam::sqr(etaP_), a_), b_)
        );

    // lambda effective
    lambdaEff_ = lambda_/(scalar(1)  + epsilon_*lambda_*tr(tau_)/etaP_);

    // Stress transport equation
    fvSymmTensorMatrix tauEqn
    (
        fvm::ddt(tau_)
      + fvm::div(phi(), tau_)
     ==
        etaPEff_/lambdaEff_*twoD
      + twoSymm(C)
      - zeta_*symm(tau_ & twoD)
      - fvm::Sp(epsilon_/etaPEff_*tr(tau_) + 1/lambdaEff_, tau_)
    );

    tauEqn.relax();
    tauEqn.solve();
}


bool Foam::viscoelasticLaws::Feta_PTT::read
(
    const dictionary& dict
)
{
    viscoelasticLaw::read(dict);
    
    Feta_PTTCoeffs_ = dict.optionalSubDict(typeName + "Coeffs");
    
    Feta_PTTCoeffs_.readEntry("rho", rho_);
    Feta_PTTCoeffs_.readEntry("etaS", etaS_);
    Feta_PTTCoeffs_.readEntry("etaP", etaP_);
    Feta_PTTCoeffs_.readEntry("epsilon", epsilon_);
    Feta_PTTCoeffs_.readEntry("lambda", lambda_);
    Feta_PTTCoeffs_.readEntry("zeta", zeta_);
    Feta_PTTCoeffs_.readEntry("A", A_);
    Feta_PTTCoeffs_.readEntry("a", a_);
    Feta_PTTCoeffs_.readEntry("b", b_);

    return true;
}
// ************************************************************************* //
