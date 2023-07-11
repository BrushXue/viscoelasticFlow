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

#include "WhiteMetznerCarreauYasuda.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscoelasticLaws
    {
        defineTypeNameAndDebug(WhiteMetznerCarreauYasuda, 0);
        addToRunTimeSelectionTable(viscoelasticLaw, WhiteMetznerCarreauYasuda, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscoelasticLaws::WhiteMetznerCarreauYasuda::WhiteMetznerCarreauYasuda
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscoelasticLaw(name, dict, U, phi),
    
    WhiteMetznerCarreauYasudaCoeffs_
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
    rho_("rho", dimDensity, WhiteMetznerCarreauYasudaCoeffs_),
    etaS_("etaS", dimDynamicViscosity, WhiteMetznerCarreauYasudaCoeffs_),
    etaP_("etaP", dimDynamicViscosity, WhiteMetznerCarreauYasudaCoeffs_),
    lambda_("lambda", dimTime, WhiteMetznerCarreauYasudaCoeffs_),
    m_("m", dimless, WhiteMetznerCarreauYasudaCoeffs_),
    n_("n", dimless, WhiteMetznerCarreauYasudaCoeffs_),
    K_("K", dimTime, WhiteMetznerCarreauYasudaCoeffs_),
    L_("L", dimTime, WhiteMetznerCarreauYasudaCoeffs_),
    a_("a", dimless, WhiteMetznerCarreauYasudaCoeffs_),
    b_("b", dimless, WhiteMetznerCarreauYasudaCoeffs_)
    
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::viscoelasticLaws::WhiteMetznerCarreauYasuda::divTau(volVectorField& U) const
{
    return
    (
        fvc::div(tau_ / rho_, "div(tau)")
      - fvc::div(etaP_ / rho_ * fvc::grad(U), "div(grad(U))")
      + fvm::laplacian((etaP_ + etaS_) / rho_, U, "laplacian(eta,U)")
    );
}


void Foam::viscoelasticLaws::WhiteMetznerCarreauYasuda::correct()
{
    // Velocity gradient tensor
    volTensorField L(fvc::grad(U()));

    // Convected derivate term
    volTensorField C(tau_ & L);

    // Twice the rate of deformation tensor
    volSymmTensorField twoD(twoSymm(L));

    // Effective viscosity and relaxation time
    volScalarField etaPValue(etaP_*Foam::pow(1.0 + Foam::pow(K_* sqrt(2.0)*mag(symm(L)),a_), (m_- 1.0)/a_));

    volScalarField lambdaValue(lambda_*Foam::pow(1.0 + Foam::pow( L_* sqrt(2.0)*mag(symm(L)),b_), (n_- 1.0)/b_));


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

bool Foam::viscoelasticLaws::WhiteMetznerCarreauYasuda::read
(
    const dictionary& dict
)
{
    viscoelasticLaw::read(dict);

    WhiteMetznerCarreauYasudaCoeffs_ = dict.optionalSubDict(typeName + "Coeffs");
    
    WhiteMetznerCarreauYasudaCoeffs_.readEntry("rho", rho_);
    WhiteMetznerCarreauYasudaCoeffs_.readEntry("etaS", etaS_);
    WhiteMetznerCarreauYasudaCoeffs_.readEntry("etaP", etaP_);
    WhiteMetznerCarreauYasudaCoeffs_.readEntry("lambda", lambda_);
    WhiteMetznerCarreauYasudaCoeffs_.readEntry("m", m_);
    WhiteMetznerCarreauYasudaCoeffs_.readEntry("n", n_);
    WhiteMetznerCarreauYasudaCoeffs_.readEntry("K", K_);
    WhiteMetznerCarreauYasudaCoeffs_.readEntry("L", L_);
    WhiteMetznerCarreauYasudaCoeffs_.readEntry("a", a_);
    WhiteMetznerCarreauYasudaCoeffs_.readEntry("b", b_);
        
    return true;
}
// ************************************************************************* //
