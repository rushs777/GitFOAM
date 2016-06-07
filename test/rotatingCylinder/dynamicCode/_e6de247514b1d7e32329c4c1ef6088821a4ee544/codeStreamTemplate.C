/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 26 "/home/simon/OpenFOAM/simon-v3.0+/run/test/rotatingCylinder/system/blockMeshDict.#codeStream"
#include "pointField.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    void codeStream_e6de247514b1d7e32329c4c1ef6088821a4ee544
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 31 "/home/simon/OpenFOAM/simon-v3.0+/run/test/rotatingCylinder/system/blockMeshDict.#codeStream"
pointField points(9);
        points[0]  = point(0.5, -0.5, -0.5);
        points[1]  = point(-0.5, -0.5, -0.5);
        points[2]  = point(0.5, 0.5, -0.5);
        points[3]  = point(-0.5, 0.5, -0.5);
	points[4]  = point(-0.5, 0, -0.5);
	points[5]  = point(0, 0, -0.5);
	points[6]  = point(0.5, 0, -0.5);
	points[7]  = point(0, 0.5, -0.5);
	points[8]  = point(0, -0.5, -0.5);


        // Duplicate z points
        label sz = points.size();
        points.setSize(2*sz);
        for (label i = 0; i < sz; i++)
        {
            const point& pt = points[i];
            points[i+sz] = point(pt.x(), pt.y(), -pt.z());
        }

        os  << points;
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

