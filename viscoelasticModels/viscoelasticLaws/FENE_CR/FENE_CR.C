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

#include "FENE_CR.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscoelasticLaws
    {
        defineTypeNameAndDebug(FENE_CR, 0);
        addToRunTimeSelectionTable(viscoelasticLaw, FENE_CR, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscoelasticLaws::FENE_CR::FENE_CR
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscoelasticLaw(name, dict, U, phi),

    FENE_CRCoeffs_
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

    rho_("rho", dimDensity, FENE_CRCoeffs_),
    etaS_("etaS", dimDynamicViscosity, FENE_CRCoeffs_),
    etaP_("etaP", dimDynamicViscosity, FENE_CRCoeffs_),
    L2_("L2", dimless, FENE_CRCoeffs_),
    lambda_("lambda", dimTime, FENE_CRCoeffs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::viscoelasticLaws::FENE_CR::divTau(volVectorField& U) const
{
    dimensionedScalar etaPEff = etaP_;

    return
    (
        fvc::div(tau_/rho_, "div(tau)")
      - fvc::laplacian(etaPEff/rho_, U, "laplacian(etaPEff,U)")
      + fvm::laplacian( (etaPEff + etaS_)/rho_, U, "laplacian(etaPEff+etaS,U)")
    );
}


void Foam::viscoelasticLaws::FENE_CR::correct()
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
        ((L2_ / lambda_ + tr(tau_)/etaP_)/(L2_ - 3.0))*etaP_*twoD
      + twoSymm(C)
      - fvm::Sp((L2_/lambda_ + tr(tau_)/etaP_)/(L2_ - 3), tau_)
    );

    tauEqn.relax();
    tauEqn.solve();
}


bool Foam::viscoelasticLaws::FENE_CR::read
(
    const dictionary& dict
)
{
    viscoelasticLaw::read(dict);
    
    FENE_CRCoeffs_ = dict.optionalSubDict(typeName + "Coeffs");
    
    FENE_CRCoeffs_.readEntry("rho", rho_);
    FENE_CRCoeffs_.readEntry("etaS", etaS_);
    FENE_CRCoeffs_.readEntry("etaP", etaP_);
    FENE_CRCoeffs_.readEntry("L2", L2_);
    FENE_CRCoeffs_.readEntry("lambda", lambda_);
    
    return true;
}
// ************************************************************************* //
