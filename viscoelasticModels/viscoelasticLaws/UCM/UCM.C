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

#include "UCM.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscoelasticLaws
    {
        defineTypeNameAndDebug(UCM, 0);
        addToRunTimeSelectionTable(viscoelasticLaw, UCM, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscoelasticLaws::UCM::UCM
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscoelasticLaw(name, dict, U, phi),

    UCMCoeffs_
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
    
    rho_("rho", dimDensity, UCMCoeffs_),
    etaP_("etaP", dimDynamicViscosity, UCMCoeffs_),
    lambda_("lambda", dimTime, UCMCoeffs_),

    etaStab_
    (
        dimensionedScalar
        (
            "etaStab",
            dimDynamicViscosity,
            dict.lookupOrDefault<scalar>("etaStab", 0.0)
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::viscoelasticLaws::UCM::divTau(volVectorField& U) const
{
    if(etaStab_.value() < SMALL)
    {
        return
        (
            fvc::div(tau_ / rho_, "div(tau)")
          - fvc::div(etaP_ / rho_ * fvc::grad(U), "div(grad(U))")
          + fvm::laplacian(etaP_ / rho_, U, "laplacian(eta,U)")
        );
    }
    else
    {
        return
        (
            fvc::div(tau_ / rho_, "div(tau)")
          - fvc::div(etaStab_ / rho_ * fvc::grad(U), "div(grad(U))")
          + fvm::laplacian(etaStab_ / rho_, U, "laplacian(eta,U)")
        );
    }
}


void Foam::viscoelasticLaws::UCM::correct()
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
      - fvm::Sp(1.0/lambda_, tau_)
    );

    tauEqn.relax();
    tauEqn.solve();
}

bool Foam::viscoelasticLaws::UCM::read
(
    const dictionary& dict
)
{
    viscoelasticLaw::read(dict);
    
    UCMCoeffs_ = dict.optionalSubDict(typeName + "Coeffs");
    
    UCMCoeffs_.readEntry("rho", rho_);
    UCMCoeffs_.readEntry("etaP", etaP_);
    UCMCoeffs_.readEntry("lambda", lambda_);
    
    return true;
}
// ************************************************************************* //
