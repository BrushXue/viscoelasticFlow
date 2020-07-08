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

#include "multiMode.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscoelasticLaws
    {
        defineTypeNameAndDebug(multiMode, 0);
        addToRunTimeSelectionTable(viscoelasticLaw, multiMode, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscoelasticLaws::multiMode::multiMode
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscoelasticLaw(name, dict, U, phi),

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

    models_()
    {
        PtrList<entry> modelEntries(dict.lookup("models"));
        models_.setSize(modelEntries.size());

        forAll (models_, modelI)
        {
            models_.set
            (
                modelI,
                viscoelasticLaw::New
                (
                    modelEntries[modelI].keyword(),
                    modelEntries[modelI].dict(),
                    U,
                    phi
                )
            );
        }
    }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::viscoelasticLaws::multiMode::divTau(volVectorField& U) const
{
    tmp<fvVectorMatrix> divMatrix = models_[0].divTau(U);
    
    for (label i = 1; i < models_.size(); i++)
    {
        divMatrix.ref() += models_[i].divTau(U);
    }

    return divMatrix;
}


Foam::tmp<Foam::volSymmTensorField> Foam::viscoelasticLaws::multiMode::tau() const
{
    tau_ *= 0.0;

    for (label i = 0; i < models_.size(); i++)
    {
        tau_ += models_[i].tau();
    }

    return tau_;
}


void Foam::viscoelasticLaws::multiMode::correct()
{
    forAll (models_, i)
    {
        Info<< "Model mode "  << i+1 << endl;
        models_[i].correct();
    }

    tau();
}

bool Foam::viscoelasticLaws::multiMode::read
(
    const dictionary& dict
)
{
    viscoelasticLaw::read(dict);
    
    return true;
}
// ************************************************************************* //
