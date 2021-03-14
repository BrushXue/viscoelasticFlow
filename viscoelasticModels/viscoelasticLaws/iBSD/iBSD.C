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

#include "iBSD.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscoelasticLaws
    {
        defineTypeNameAndDebug(iBSD, 0);
        addToRunTimeSelectionTable(viscoelasticLaw, iBSD, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscoelasticLaws::iBSD::iBSD
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscoelasticLaw(name, dict, U, phi),

    iBSDCoeffs_
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

    rho_("rho", dimDensity, iBSDCoeffs_),
    etaS_("etaS", dimDynamicViscosity, iBSDCoeffs_),
    etaP_("etaP", dimDynamicViscosity, iBSDCoeffs_),
    lambda_("lambda", dimTime, iBSDCoeffs_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::viscoelasticLaws::iBSD::divTau(volVectorField& U) const
{
    volTensorField gradU(fvc::grad(etaP_ * U, "grad(etaP,U)"));
    return
    (
        fvc::div(tau_/rho_, "div(tau)")
      - fvc::div(gradU/rho_, "div(gradU)")
      + fvm::laplacian( (etaP_ + etaS_)/rho_, U, "laplacian(etaP+etaS,U)")
    );
}


void Foam::viscoelasticLaws::iBSD::correct()
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
      - fvm::Sp(1/lambda_, tau_)
    );

    tauEqn.relax();
    tauEqn.solve();
}


bool Foam::viscoelasticLaws::iBSD::read
(
    const dictionary& dict
)
{
    viscoelasticLaw::read(dict);
    
    iBSDCoeffs_ = dict.optionalSubDict(typeName + "Coeffs");
    
    iBSDCoeffs_.readEntry("rho", rho_);
    iBSDCoeffs_.readEntry("etaS", etaS_);
    iBSDCoeffs_.readEntry("etaP", etaP_);
    iBSDCoeffs_.readEntry("lambda", lambda_);
    
    return true;
}
// ************************************************************************* //
