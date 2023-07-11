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

#include "FENE_CR_LCR.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscoelasticLaws
    {
        defineTypeNameAndDebug(FENE_CR_LCR, 0);
        addToRunTimeSelectionTable(viscoelasticLaw, FENE_CR_LCR, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscoelasticLaws::FENE_CR_LCR::FENE_CR_LCR
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscoelasticLaw(name, dict, U, phi),

    FENE_CR_LCR_Coeffs_
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
    
    Psi
    (
        IOobject
        (
            "Psi" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    E
    (
        IOobject
        (
            "E" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedSymmTensor
        (
                "I",
                dimless,
                pTraits<symmTensor>::I
        ),
         extrapolatedCalculatedFvPatchField<symmTensor>::typeName
    ),

    R
    (
        IOobject
        (
            "R" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedTensor
        (
                "I",
                dimless,
                pTraits<tensor>::I
        ),
         extrapolatedCalculatedFvPatchField<tensor>::typeName
    ),

    rho_("rho", dimDensity, FENE_CR_LCR_Coeffs_),
    etaS_("etaS", dimDynamicViscosity, FENE_CR_LCR_Coeffs_),
    etaP_("etaP", dimDynamicViscosity, FENE_CR_LCR_Coeffs_),
    L2_("L2", dimless, FENE_CR_LCR_Coeffs_),
    lambda_("lambda", dimTime, FENE_CR_LCR_Coeffs_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::viscoelasticLaws::FENE_CR_LCR::divTau(volVectorField& U) const
{
    return
    (
        fvc::div(tau_ / rho_, "div(tau)")
      - fvc::div(etaP_ / rho_ * fvc::grad(U), "div(grad(U))")
      + fvm::laplacian((etaP_ + etaS_) / rho_, U, "laplacian(eta,U)")
    );
}

void Foam::viscoelasticLaws::FENE_CR_LCR::correct()
{
    // OpenFOAM grad(U) is transposed
    const volTensorField L(fvc::grad(U()));
    
    // Rotate gradU into the principal axes
    const volTensorField M(R.T() & L.T() & R);
    
    // Calculate S
    volSymmTensorField S(twoSymm(M));
    decompose(M, E, R, S);
    
    // Calculate extensibility
    volScalarField f(L2_ / (L2_ - tr(R & E & R.T())));
    
    // Calculate F
    volSymmTensorField F((inv(E) - symmTensor::I) * f / lambda_);

    // Solve Psi
    fvSymmTensorMatrix PsiEqn
    (
        fvm::ddt(Psi)
      + fvm::div(phi(), Psi)
     ==
        symm(R & (S + F) & R.T())
    );
    PsiEqn.relax();
    PsiEqn.solve();
    
    // Calculate eigenvalues
    eigen(Psi, E, R);

    // Reconstruct stress
    tau_=etaP_ / lambda_ * f * (symm(R & E & R.T()) - symmTensor::I);
    tau_.correctBoundaryConditions();
}

bool Foam::viscoelasticLaws::FENE_CR_LCR::read
(
    const dictionary& dict
)
{
    viscoelasticLaw::read(dict);
    
    FENE_CR_LCR_Coeffs_ = dict.optionalSubDict(typeName + "Coeffs");
    
    FENE_CR_LCR_Coeffs_.readEntry("rho", rho_);
    FENE_CR_LCR_Coeffs_.readEntry("etaS", etaS_);
    FENE_CR_LCR_Coeffs_.readEntry("etaP", etaP_);
    FENE_CR_LCR_Coeffs_.readEntry("L2", L2_);
    FENE_CR_LCR_Coeffs_.readEntry("lambda", lambda_);
    
    return true;
}
// ************************************************************************* //
