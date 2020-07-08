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

#include "S_MDCPP.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscoelasticLaws
    {
        defineTypeNameAndDebug(S_MDCPP, 0);
        addToRunTimeSelectionTable(viscoelasticLaw, S_MDCPP, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscoelasticLaws::S_MDCPP::S_MDCPP
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscoelasticLaw(name, dict, U, phi),

    S_MDCPPCoeffs_
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
    I_
    (
        dimensionedSymmTensor
        (
            "I",
            dimless,
            symmTensor
            (
                1, 0, 0,
                   1, 0,
                      1
            )
        )
    ),
    rho_("rho", dimDensity, S_MDCPPCoeffs_),
    etaS_("etaS", dimDynamicViscosity, S_MDCPPCoeffs_),
    etaP_("etaP", dimDynamicViscosity, S_MDCPPCoeffs_),
    zeta_("zeta", dimless, S_MDCPPCoeffs_),
    lambdaOb_("lambdaOb", dimTime, S_MDCPPCoeffs_),
    lambdaOs_("lambdaOs", dimTime, S_MDCPPCoeffs_),
    q_("q", dimless, S_MDCPPCoeffs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::viscoelasticLaws::S_MDCPP::divTau(volVectorField& U) const
{
    dimensionedScalar etaPEff = etaP_;

    return
    (
        fvc::div(tau_/rho_, "div(tau)")
      - fvc::laplacian(etaPEff/rho_, U, "laplacian(etaPEff,U)")
      + fvm::laplacian( (etaPEff + etaS_)/rho_, U, "laplacian(etaPEff+etaS,U)")
    );

}


void Foam::viscoelasticLaws::S_MDCPP::correct()
{
    // Velocity gradient tensor
    volTensorField L(fvc::grad(U()));

    // Convected derivate term
    volTensorField C(tau_ & L);

    // Twice the rate of deformation tensor
    volSymmTensorField twoD(twoSymm(L));

    // Lambda (Backbone stretch)
    volScalarField Lambda(Foam::sqrt(1 + tr(tau_)*lambdaOb_*(1 - zeta_)/3/etaP_));

    // Auxiliary field
    volScalarField aux(Foam::exp( 2/q_*(Lambda - 1)));

    // Extra function
    volScalarField fTau(aux*(2*lambdaOb_/lambdaOs_*(1 - 1/Lambda) + 1/Foam::sqr(Lambda)));

    // Stress transport equation
    fvSymmTensorMatrix tauEqn
    (
        fvm::ddt(tau_)
      + fvm::div(phi(), tau_)
     ==
        etaP_/lambdaOb_*twoD
      + twoSymm(C)
      - zeta_*symm(tau_ & twoD)
      - fvm::Sp(1/lambdaOb_*fTau, tau_)
      - (
            1/lambdaOb_*(etaP_/lambdaOb_/(1 - zeta_)*(fTau - aux)*I_)
        )
    );

    tauEqn.relax();
    tauEqn.solve();
}

bool Foam::viscoelasticLaws::S_MDCPP::read
(
    const dictionary& dict
)
{
    viscoelasticLaw::read(dict);
    
    S_MDCPPCoeffs_ = dict.optionalSubDict(typeName + "Coeffs");
    
    S_MDCPPCoeffs_.readEntry("rho", rho_);
    S_MDCPPCoeffs_.readEntry("etaS", etaS_);
    S_MDCPPCoeffs_.readEntry("etaP", etaP_);
    S_MDCPPCoeffs_.readEntry("zeta", zeta_);
    S_MDCPPCoeffs_.readEntry("lambdaOb", lambdaOb_);
    S_MDCPPCoeffs_.readEntry("lambdaOs", lambdaOs_);
    S_MDCPPCoeffs_.readEntry("q", q_);
    
    return true;
}// ************************************************************************* //
