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

#include "DCPP.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscoelasticLaws
    {
        defineTypeNameAndDebug(DCPP, 0);
        addToRunTimeSelectionTable(viscoelasticLaw, DCPP, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscoelasticLaws::DCPP::DCPP
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscoelasticLaw(name, dict, U, phi),

    DCPPCoeffs_
    (
        dict.optionalSubDict(typeName + "Coeffs")
    ),

    S_
    (
        IOobject
        (
            "S" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    Lambda_
    (
        IOobject
        (
            "Lambda" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    tau_
    (
        IOobject
        (
            "tau" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedSymmTensor
        (
            "zero",
            dimPressure,
            symmTensor::zero
        )
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

    rho_("rho", dimDensity, DCPPCoeffs_),
    etaS_("etaS", dimDynamicViscosity, DCPPCoeffs_),
    etaP_("etaP", dimDynamicViscosity, DCPPCoeffs_),
    zeta_("zeta", dimless, DCPPCoeffs_),
    lambdaOb_("lambdaOb", dimTime, DCPPCoeffs_),
    lambdaOs_("lambdaOs", dimTime, DCPPCoeffs_),
    q_("q", dimless, DCPPCoeffs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::viscoelasticLaws::DCPP::divTau(volVectorField& U) const
{
    dimensionedScalar etaPEff = etaP_;

    return
    (
        fvc::div(tau_/rho_, "div(tau)")
      - fvc::laplacian(etaPEff/rho_, U, "laplacian(etaPEff,U)")
      + fvm::laplacian((etaPEff + etaS_)/rho_, U, "laplacian(etaPEff+etaS,U)")
    );
}


void Foam::viscoelasticLaws::DCPP::correct()
{
    // Velocity gradient tensor
    volTensorField L(fvc::grad(U()));

    // Upper convected derivate term
    volTensorField Cupper(S_ & L);

    // Lower convected derivate term
    volTensorField Clower(L & S_);

    // Twice the rate of deformation tensor
    volSymmTensorField twoD(twoSymm(L));

    // Evolution of orientation
    fvSymmTensorMatrix SEqn
    (
        fvm::ddt(S_)
      + fvm::div(phi(), S_)
     ==
        (1 - zeta_/2)*twoSymm(Cupper)
      - (zeta_/2)*twoSymm(Clower)
      - (1 - zeta_)*fvm::Sp((twoD && S_), S_)
      - fvm::Sp(1/lambdaOb_/Foam::sqr(Lambda_), S_)
      + 1/lambdaOb_/Foam::sqr(Lambda_)/3*I_
    );

    SEqn.relax();
    SEqn.solve();

     // Evolution of the backbone stretch
    fvScalarMatrix lambdaEqn
    (
        fvm::ddt(Lambda_)
      + fvm::div(phi(), Lambda_)
     ==
        fvm::Sp((twoD && S_)/2 , Lambda_)
      - fvm::Sp(Foam::exp(2/q_*(Lambda_ - 1))/lambdaOs_ , Lambda_)
      + Foam::exp( 2/q_*(Lambda_ - 1))/lambdaOs_
    );

    lambdaEqn.relax();
    lambdaEqn.solve();

    // Viscoelastic stress
    tau_ = etaP_/lambdaOb_/(1 - zeta_) * (3*Foam::sqr(Lambda_)*S_ - I_);
}


bool Foam::viscoelasticLaws::DCPP::read
(
    const dictionary& dict
)
{
    viscoelasticLaw::read(dict);
    
    DCPPCoeffs_ = dict.optionalSubDict(typeName + "Coeffs");
    
    DCPPCoeffs_.readEntry("rho", rho_);
    DCPPCoeffs_.readEntry("etaS", etaS_);
    DCPPCoeffs_.readEntry("etaP", etaP_);
    DCPPCoeffs_.readEntry("zeta", zeta_);
    DCPPCoeffs_.readEntry("lambdaOb", lambdaOb_);
    DCPPCoeffs_.readEntry("lambdaOs", lambdaOs_);
    DCPPCoeffs_.readEntry("q", q_);
    
    return true;
}
// ************************************************************************* //
