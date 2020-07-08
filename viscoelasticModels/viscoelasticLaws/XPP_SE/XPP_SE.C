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

#include "XPP_SE.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscoelasticLaws
    {
        defineTypeNameAndDebug(XPP_SE, 0);
        addToRunTimeSelectionTable(viscoelasticLaw, XPP_SE, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscoelasticLaws::XPP_SE::XPP_SE
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscoelasticLaw(name, dict, U, phi),

    XPP_SECoeffs_
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
    rho_("rho", dimDensity, XPP_SECoeffs_),
    etaS_("etaS", dimDynamicViscosity, XPP_SECoeffs_),
    etaP_("etaP", dimDynamicViscosity, XPP_SECoeffs_),
    alpha_("alpha", dimless, XPP_SECoeffs_),
    lambdaOb_("lambdaOb", dimTime, XPP_SECoeffs_),
    lambdaOs_("lambdaOs", dimTime, XPP_SECoeffs_),
    q_("q", dimless, XPP_SECoeffs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::viscoelasticLaws::XPP_SE::divTau(volVectorField& U) const
{

    dimensionedScalar etaPEff = etaP_;

    return
    (
        fvc::div(tau_/rho_, "div(tau)")
      - fvc::laplacian(etaPEff/rho_, U, "laplacian(etaPEff,U)")
      + fvm::laplacian( (etaPEff + etaS_)/rho_, U, "laplacian(etaPEff+etaS,U)")
    );

}


void Foam::viscoelasticLaws::XPP_SE::correct()
{
    // Velocity gradient tensor
    volTensorField L(fvc::grad(U()));

    // Convected derivate term
    volTensorField C(tau_ & L);

    // Twice the rate of deformation tensor
    volSymmTensorField twoD(twoSymm(L));

    // Lambda (Backbone stretch)
    volScalarField Lambda(Foam::sqrt(1 + tr(tau_)*lambdaOb_/3/etaP_));

    // lambdaS (stretch relaxation time)
    volScalarField lambdaS(lambdaOs_*Foam::exp(-2/q_*(Lambda - 1)));

    // Extra function
    volScalarField fTau(2*lambdaOb_/lambdaS*(1 - 1/Lambda)
      + 1/Foam::sqr(Lambda)*
        (1 - alpha_*tr(tau_ & tau_)/3/Foam::sqr(etaP_/lambdaOb_)));

     // Stress transport equation
    fvSymmTensorMatrix tauEqn
    (
        fvm::ddt(tau_)
      + fvm::div(phi(), tau_)
     ==
        etaP_/lambdaOb_*twoD
      + twoSymm(C)
      - fvm::Sp(1/lambdaOb_*fTau, tau_)
      - (
            1/lambdaOb_*
            (
                alpha_*lambdaOb_/etaP_*symm(tau_ & tau_)
              + etaP_/lambdaOb_*(fTau - 1)*I_
            )
        )
    );

    tauEqn.relax();
    tauEqn.solve();

}


bool Foam::viscoelasticLaws::XPP_SE::read
(
    const dictionary& dict
)
{
    viscoelasticLaw::read(dict);
    
    XPP_SECoeffs_ = dict.optionalSubDict(typeName + "Coeffs");
    
    XPP_SECoeffs_.readEntry("rho", rho_);
    XPP_SECoeffs_.readEntry("etaS", etaS_);
    XPP_SECoeffs_.readEntry("etaP", etaP_);
    XPP_SECoeffs_.readEntry("alpha", alpha_);
    XPP_SECoeffs_.readEntry("lambdaOb", lambdaOb_);
    XPP_SECoeffs_.readEntry("lambdaOs", lambdaOs_);
    XPP_SECoeffs_.readEntry("q", q_);
    
    return true;
}
// ************************************************************************* //
