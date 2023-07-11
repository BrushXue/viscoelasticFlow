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

#include "WhiteMetznerCross.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscoelasticLaws
    {
        defineTypeNameAndDebug(WhiteMetznerCross, 0);
        addToRunTimeSelectionTable(viscoelasticLaw, WhiteMetznerCross, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscoelasticLaws::WhiteMetznerCross::WhiteMetznerCross
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscoelasticLaw(name, dict, U, phi),

    WhiteMetznerCrossCoeffs_
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
    rho_("rho", dimDensity, WhiteMetznerCrossCoeffs_),
    etaS_("etaS", dimDynamicViscosity, WhiteMetznerCrossCoeffs_),
    etaP_("etaP", dimDynamicViscosity, WhiteMetznerCrossCoeffs_),
    lambda_("lambda", dimTime, WhiteMetznerCrossCoeffs_),
    K_("K", dimTime, WhiteMetznerCrossCoeffs_),
    L_("L", dimTime, WhiteMetznerCrossCoeffs_),
    m_("m", dimless, WhiteMetznerCrossCoeffs_),
    n_("n", dimless, WhiteMetznerCrossCoeffs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::viscoelasticLaws::WhiteMetznerCross::divTau(volVectorField& U) const
{
    return
    (
        fvc::div(tau_ / rho_, "div(tau)")
      - fvc::div(etaP_ / rho_ * fvc::grad(U), "div(grad(U))")
      + fvm::laplacian((etaP_ + etaS_) / rho_, U, "laplacian(eta,U)")
    );
}


void Foam::viscoelasticLaws::WhiteMetznerCross::correct()
{
    // Velocity gradient tensor
    volTensorField L(fvc::grad(U()));

    // Convected derivate term
    volTensorField C(tau_ & L);

    // Twice the rate of deformation tensor
    volSymmTensorField twoD(twoSymm(L));

    // Effective viscosity and relaxation time
    volScalarField etaPValue(etaP_/(1.0 + Foam::pow(K_* sqrt(2.0)*mag(symm(L)), (1.0 - m_))));

    volScalarField lambdaValue(lambda_/(1.0 + Foam::pow(L_ * sqrt(2.0)*mag(symm(L)), (1.0 - n_))));

    // Stress transport equation
    fvSymmTensorMatrix tauEqn
    (
        fvm::ddt(tau_)
      + fvm::div(phi(), tau_)
     ==
        etaPValue/lambdaValue*twoD
      + twoSymm(C)
      - fvm::Sp(1.0/lambdaValue, tau_)
    );

    tauEqn.relax();
    tauEqn.solve();
}


bool Foam::viscoelasticLaws::WhiteMetznerCross::read
(
    const dictionary& dict
)
{
    viscoelasticLaw::read(dict);

    WhiteMetznerCrossCoeffs_ = dict.optionalSubDict(typeName + "Coeffs");
    
    WhiteMetznerCrossCoeffs_.readEntry("rho", rho_);
    WhiteMetznerCrossCoeffs_.readEntry("etaS", etaS_);
    WhiteMetznerCrossCoeffs_.readEntry("etaP", etaP_);
    WhiteMetznerCrossCoeffs_.readEntry("lambda", lambda_);
    WhiteMetznerCrossCoeffs_.readEntry("K", K_);
    WhiteMetznerCrossCoeffs_.readEntry("L", L_);
    WhiteMetznerCrossCoeffs_.readEntry("m", m_);
    WhiteMetznerCrossCoeffs_.readEntry("n", n_);

    return true;
}

// ************************************************************************* //
