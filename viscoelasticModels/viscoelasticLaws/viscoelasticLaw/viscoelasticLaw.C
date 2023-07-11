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

#include "viscoelasticLaw.H"
#include <Eigenvalues>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(viscoelasticLaw, 0);
    defineRunTimeSelectionTable(viscoelasticLaw, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

viscoelasticLaw::viscoelasticLaw
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    name_(name),
    dict_(dict),
    U_(U),
    phi_(phi)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool viscoelasticLaw::read(const dictionary& dict)
{
    dict_ = dict;
    
    return true;
}

void viscoelasticLaw::decompose
(
    const volTensorField& M,
    const volSymmTensorField& E, 
    const volTensorField& R,
    volSymmTensorField& S
)
{
    const scalar epsilon=1e-9;
    forAll(S, celli) {
        symmTensor& tS = S[celli];
        const symmTensor& tE = E[celli];
        const tensor& tM = M[celli];
        if (mag(tE.xx() - tE.yy()) > epsilon)
            tS.xy() = (tE.yy() * tM.xy() + tE.xx() * tM.yx()) * 
            (log(tE.yy()) - log(tE.xx())) / (tE.yy() - tE.xx());
        if (mag(tE.xx() - tE.zz()) > epsilon)
            tS.xz() = (tE.zz() * tM.xz() + tE.xx() * tM.zx()) * 
            (log(tE.zz()) - log(tE.xx())) / (tE.zz() - tE.xx());
        if (mag(tE.yy() - tE.zz()) > epsilon)
            tS.yz() = (tE.zz() * tM.yz() + tE.yy() * tM.zy()) * 
            (log(tE.zz()) - log(tE.yy())) / (tE.zz() - tE.yy());
    }
}

void viscoelasticLaw::eigen
(
    const volSymmTensorField& Psi,
    volSymmTensorField& E,
    volTensorField& R
)
{
    Eigen::Matrix3d A;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> sol;
    // Update eigen (dimless)
    forAll(Psi, celli) {
        const symmTensor& tPsi = Psi[celli]; 
        A(0,0)=tPsi.xx();
        A(1,1)=tPsi.yy();
        A(2,2)=tPsi.zz();
        A(1,0)=tPsi.xy();
        A(2,0)=tPsi.xz();
        A(2,1)=tPsi.yz();
        sol.compute(A);
        Eigen::Vector3d val = sol.eigenvalues();
        Eigen::Matrix3d vec = sol.eigenvectors();
        symmTensor& tE = E[celli];
        tE = symmTensor(exp(val(0)),0,0,exp(val(1)),0,exp(val(2)));
        /*
        tE.xx() = exp(val(0));
        tE.yy() = exp(val(1));
        tE.zz() = exp(val(2));
        */
        tensor& tR = R[celli];
        tR = tensor(vec(0, 0), vec(0, 1), vec(0, 2), vec(1, 0), vec(01, 1), vec(1, 2), vec(2, 0), vec(2, 1), vec(2, 2));
        /*
        tR.xx() = vec(0, 0);
        tR.xy() = vec(0, 1);
        tR.xz() = vec(0, 2);
        tR.yx() = vec(1, 0);
        tR.yy() = vec(1, 1);
        tR.yz() = vec(1, 2);
        tR.zx() = vec(2, 0);
        tR.zy() = vec(2, 1);
        tR.zz() = vec(2, 2);
        */
    }
}

} // End namespace Foam
// ************************************************************************* //
