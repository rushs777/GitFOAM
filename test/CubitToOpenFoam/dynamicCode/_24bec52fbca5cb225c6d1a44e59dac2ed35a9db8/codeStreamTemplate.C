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
    void codeStream_24bec52fbca5cb225c6d1a44e59dac2ed35a9db8
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 31 "/home/simon/OpenFOAM/simon-v3.0+/run/test/rotatingCylinder/system/blockMeshDict.#codeStream"
pointField points(32);
	// Southwest points
	points[0]  = point(-1.0, -1.0, -1.0);
	points[1]  = point(-0.25, -1.0, -1.0);
	points[2]  = point(0.0, -1.0, -1.0);
	points[3]  = point(0.0, -0.5, -1.0);
	points[4]  = point(0.0, -0.25, -1.0);
	points[5]  = point(-0.125, -0.2165064, -1.0);
	points[6]  = point(-0.25, 0.0, -1.0);
	points[7]  = point(-0.5, 0.0, -1.0);
	points[8]  = point(-1.0, 0.0, -1.0);
	points[9]  = point(-1.0, -0.4330127, -1.0);
	points[10] = point(-0.25, -0.4330127, -1.0);
	// Southeast points
	points[11] = point(0.25, -1.0, -1.0);
	points[12] = point(1.0, -1.0, -1.0);
	points[13] = point(1.0, -0.4330127, -1.0);
	points[14] = point(1.0, 0.0, -1.0);
	points[15] = point(0.5, 0.0, -1.0);
	points[16] = point(0.25, 0.0, -1.0);
	points[17] = point(0.125, -0.2165064, -1.0);
	points[18] = point(0.25, -0.4330127, -1.0);
	// Northeast points
	points[19] = point(1.0, 0.4330127, -1.0);
	points[20] = point(1.0, 1.0, -1.0);
	points[21] = point(0.4330127, 1.0, -1.0);
	points[22] = point(0.0, 1.0, -1.0);
	points[23] = point(0.0, 0.5, -1.0);
	points[24] = point(0.0, 0.25, -1.0);
	points[25] = point(0.125, 0.2165064, -1.0);
	points[26] = point(0.25, 0.4330127, -1.0);
	// Northwest points
	points[27] = point(-0.25, 1.0, -1.0);
	points[28] = point(-1.0, 1.0, -1.0);
	points[29] = point(-1.0, 0.4330127, -1.0);
	points[30] = point(-0.125, 0.2165064, -1.0);
	points[31] = point(-0.25, 0.4330127, -1.0);


/*
        points[0]  = point(0.5, -0.5, -0.5);
        points[1]  = point(-0.5, -0.5, -0.5);
        points[2]  = point(0.5, 0.5, -0.5);
        points[3]  = point(-0.5, 0.5, -0.5);
	points[4]  = point(-0.5, 0, -0.5);
	points[5]  = point(0.5, 0, -0.5);
	points[6]  = point(0, 0.5, -0.5);
	points[7]  = point(0, -0.5, -0.5);
	// Inner Cylinder
	points[8]  = point(0, -0.1, -0.5); //South
	points[9]  = point(-0.1, 0, -0.5); // West
	points[10] = point(0, 0.1, -0.5); // North
	points[11] = point[0.1, 0, -0.5); // East
	// Outer Cylinder
	points[12]  = point(0, -0.3, -0.5); //South
	points[13]  = point(-0.3, 0, -0.5); // West
	points[14] = point(0, 0.3, -0.5); // North
	points[15] = point[0.3, 0, -0.5); // East
*/


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

