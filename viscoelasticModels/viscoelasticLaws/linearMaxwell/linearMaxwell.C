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

#include "linearMaxwell.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscoelasticLaws
    {
        defineTypeNameAndDebug(linearMaxwell, 0);
        addToRunTimeSelectionTable(viscoelasticLaw, linearMaxwell, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscoelasticLaws::linearMaxwell::linearMaxwell
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscoelasticLaw(name, dict, U, phi),

    linearMaxwellCoeffs_
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

    rho_("rho", dimDensity, linearMaxwellCoeffs_),
    etaS_("etaS", dimDynamicViscosity, linearMaxwellCoeffs_),
    etaP_("etaP", dimDynamicViscosity, linearMaxwellCoeffs_),
    lambda_("lambda", dimTime, linearMaxwellCoeffs_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::viscoelasticLaws::linearMaxwell::divTau(volVectorField& U) const
{
    dimensionedScalar etaPEff = etaP_;

    return
    (
        fvc::div(tau_/rho_, "div(tau)")
      - fvc::laplacian(etaPEff/rho_, U, "laplacian(etaPEff,U)")
      + fvm::laplacian( (etaPEff + etaS_)/rho_, U, "laplacian(etaPEff+etaS,U)")
    );
}


void Foam::viscoelasticLaws::linearMaxwell::correct()
{
    // Velocity gradient tensor
    volTensorField L(fvc::grad(U()));

    // Twice the rate of deformation tensor
    volSymmTensorField twoD(twoSymm(L));

     // Stress transport equation
    fvSymmTensorMatrix tauEqn
    (
        fvm::ddt(tau_)
     ==
        etaP_/lambda_*twoD
      - fvm::Sp( 1/lambda_, tau_)
    );

    tauEqn.relax();
    tauEqn.solve();
}

bool Foam::viscoelasticLaws::linearMaxwell::read
(
    const dictionary& dict
)
{
    viscoelasticLaw::read(dict);
    
    linearMaxwellCoeffs_ = dict.optionalSubDict(typeName + "Coeffs");
    
    linearMaxwellCoeffs_.readEntry("rho", rho_);
    linearMaxwellCoeffs_.readEntry("etaS", etaS_);
    linearMaxwellCoeffs_.readEntry("etaP", etaP_);
    linearMaxwellCoeffs_.readEntry("lambda", lambda_);
    
    return true;
}
// ************************************************************************* //
