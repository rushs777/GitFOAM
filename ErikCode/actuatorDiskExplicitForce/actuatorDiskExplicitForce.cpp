/*-------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-----------------------------------------------------------------------
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


Application
    actuatorDiskExplicitForce

Description
	Add volume force in an actuator disk region from thrust , 		torque and geometry defined in fvSolution.
	Use with actuatorDiskExplicitForceSimpleFoam

	The actuator disk can be defined by adding the following lines in fvSolution:

	actuatorDisk
	{
		interiorRadius 1.6; // Radius of the propeller hub
		exteriorRadius 20.5; // Exterior radius of the propeller
		thrust 47.5e3; // Total force in the axial direction [N]
		torque 112.0 e3; // Total torque in the actuator disk region , positive according to the right hand rule
	density 1.2; // Fluid desity
	startPoint (103.0 0 0); // Coordinates of start point
	endPoint (102.0 0 0); // Coordinates of end point
	}

	Written by Erik Svenning , October 2010


\*-------------------------------------------------------------------*/

#include "actuatorDiskExplicitForce.h"

#include "faceAreaPairGAMGAgglomeration.H"
#include "fvMesh .H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam{
defineTypeNameAndDebug(actuatorDiskExplicitForce, 0);


//Default constructor
actuatorDiskExplicitForce::actuatorDiskExplicitForce(){

	//Set default values to all member variables
	mPointStartCenterLine.x() = 0.0;
	mPointStartCenterLine.y() = 0.0;
	mPointStartCenterLine.z() = 0.0;

	mPointEndCenterLine.x() = 0.0;
	mPointEndCenterLine.y() = 0.0;
	mPointEndCenterLine.z() = 0.0;

	mExtRadius = 0.0;
	mIntRadius = 0.0;
	mThrust = 0.0;
	mTorque = 0.0;
	mRho = 1.0;
}

actuatorDiskExplicitForce::~actuatorDiskExplicitForce(){

}

void actuatorDiskExplicitForce::ReadGeometry(const fvMesh &iMesh) {
	if(debug >= 1) {
		Info << "Reading actuator disk geometry.\n";
	}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Read actuator dict definition from solution dict ( fvSolution )
//////////////////////////////////////////////////////////////////////////////////////////////////////
	Istream& is1 = iMesh.solutionDict().subDict("actuatorDisk").lookup("interiorRadius");
	is1.format(IOstream::ASCII);
	is1 >> mIntRadius;

	Istream& is2 = iMesh.solutionDict().subDict("actuatorDisk").lookup("exteriorRadius");
	is2.format(IOstream::ASCII);
	is2 >> mExtRadius;

	Istream& is3 = iMesh.solutionDict().subDict("actuatorDisk").lookup("thrust");
	is3.format(IOstream::ASCII);
	is3 >> mThrust;

	Istream& is4 = iMesh.solutionDict().subDict("actuatorDisk").lookup("torque");
	is4.format(IOstream::ASCII);
	is4 >> mTorque;

	Istream& is6 = iMesh.solutionDict().subDict("actuatorDisk").lookup("density");
	is6.format(IOstream::ASCII);
	is6 >> mRho;

	Istream& is7 = iMesh.solutionDict().subDict("actuatorDisk").lookup("startPoint");
	is7.format(IOstream::ASCII);
	is7 >> mPointStartCenterLine;

	Istream& is8 = iMesh.solutionDict().subDict("actuatorDisk").lookup("endPoint");
	is8.format(IOstream::ASCII);
	is8 >> mPointEndCenterLine;

	if(debug >= 2) {
		Info << "Actuator disk values loaded from fvSolution:\n";
		Info << "mIntRadius: " << mIntRadius << "\n";
		Info << "mExtRadius: " << mExtRadius << "\n";
		Info << "mThrust: " << mThrust << "\n";
		Info << "mTorque: " << mTorque << "\n";
		Info << "mRho: " << mRho << "\n";
		Info << "mPointStartCenterLine: " << mPointStartCenterLine << "\n";
		Info << "mPointEndCenterLine: " << mPointEndCenterLine << "\n";
	}

}

void actuatorDiskExplicitForce::CalcActuatorDiskVolForce(const fvMesh &iMesh, vectorField &ioVolumeForce) {

	if(debug >= 1) {
		Info << "Calculating volume force from actuator disk.\n";
	}

	ReadGeometry(iMesh);

	scalar RadialDist2;
	vector LineTangent;
	vector CircumferentialDirection;

	vector TotalForce(0.0,0.0,0.0);
	scalar TotalTorque = 0.0;

	scalar DiskVolume = 0;

	// Loop over all cells and check if the cell center is in the actuator disk region
	for (label i = 0; i < iMesh.C().size(); i++) {

		if(PointIsInDisk(mPointStartCenterLine, mPointEndCenterLine, iMesh.C()[i], RadialDist2, LineTangent, CircumferentialDirection)) {

			if(debug >= 3) {
				Info << "Point: " << i << " is in the actuator disk. Coordinates: " << iMesh.C()[i] << "\n";
			}

			vector axialForce = LineTangent * CalcAxialForce( sqrt ( RadialDist2 ) , mRho ) / mRho ;
			ioVolumeForce[i] += axialForce;

			// compute the total force added to the actuator disk , this is just for control
			TotalForce += axialForce * iMesh.V()[i];

			vector circForce = CircumferentialDirection * CalcCircForce ( sqrt (RadialDist2 ) , mRho ) / mRho ;
			ioVolumeForce[i] += circForce;

			TotalTorque += ( CalcCircForce ( sqrt ( RadialDist2 ) , mRho )/ mRho )* sqrt ( RadialDist2 )*iMesh.V()[i];
			DiskVolume += iMesh.V()[i];
		}
	}

	Info << "Total axial force: " << TotalForce << "\n";
	Info << "Total torque: " << TotalTorque << "\n";
	Info << "Total disk volume: " << DiskVolume << "\n";
}

void actuatorDiskExplicitForce::WriteVTK() {
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
180 // Write the outer surface of the actuator disk to a VTK file so that it can be visualized
in Paraview .
181 //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
182 FILE * file ;
183 char fileName [100];
184
185 unsigned int NumCells = 20; // The cylindrical surface is visualized as 20
rectangular surfaces
186 unsigned int NumPoints = 40; // 40 points are required ; 20 points at each end of
the cylinder
187 unsigned int NumInts = 5* NumCells ; // Number of integers needed in the the VTK
file ; each surface has 4 corner points , so we need 4 corner indices + the
index of the surface = 5 indices per surface
188
189 vectorField points ( NumPoints , vector :: zero );
190
191 vector VecLineTangent ( mPointEndCenterLine - mPointStartCenterLine );
192 scalar LineTangentLength = sqrt ( VecLineTangent .x () * VecLineTangent .x () +
VecLineTangent .y () * VecLineTangent .y () + VecLineTangent . z () * VecLineTangent .
z () );
193
194 if( LineTangentLength != 0.0) {
195 VecLineTangent /= LineTangentLength ;
196 }
197 else {
198 Info << " Warning : The centerline tangent has zero length .\n";
199 return ;
200 }
201
202 // We need to find a vector in the radial direction . This can be any vector as
long as it points in the radial direction .
203 // First try with (1 0 0) and see if we can project it onto the normal plane of
the actuator disk resulting in a vector in
204 // the radial direction .
205 vector VecRadialDirection (1.0 ,0.0 ,0.0) ;
206 VecRadialDirection -= ( VecRadialDirection & VecLineTangent ) * VecLineTangent ;
207
208 if( mag ( VecRadialDirection ) < SMALL ) {
209 // If we enter this if statement , our guess (1 0 0) was parallel to the
centerline of the actuator disk . Then
210 // we try (0 1 0) instead . Since (1 0 0) was parallel to the centerline ,
(0 1 0) will for sure not be parallel to
211 // the centerline .
212 VecRadialDirection . x () = 0.0;
213 VecRadialDirection . y () = 1.0;
214 VecRadialDirection . z () = 0.0;
215
216 VecRadialDirection -= ( VecRadialDirection & VecLineTangent ) *
VecLineTangent ;
217 }
218
219 if( mag ( VecRadialDirection ) > SMALL ) {
220 VecRadialDirection /= mag ( VecRadialDirection );
221 }
222 else {
223 Info << " Warning in actuatorDiskExplicitForce :: WriteVTK (): mag (
VecRadialDirection ) close to zero .\n";
224 }
225
226 vector VecRadialDirection2 = VecLineTangent ^ VecRadialDirection ;
227 scalar XLocal = 0.0 , YLocal = 0.0;
228
229 // Compute points on first side of disk region
230 double phi = 0.0;
231 for ( unsigned int i = 0; i < NumCells ; i ++) {
232 XLocal = mExtRadius * cos ( phi );
233 YLocal = mExtRadius * sin ( phi );
234
235 vector point ( mPointStartCenterLine + XLocal * VecRadialDirection + YLocal
* VecRadialDirection2 ) ;
236 points [ i] = point ;
237 phi += (1.0/ double ( NumCells )) *2* mPI ;
238 }
239
240 // Compute points on second side of disk region
241 phi = 0.0;
242 for ( unsigned int i = 0; i < NumCells ; i ++) {
243 XLocal = mExtRadius * cos ( phi );
244 YLocal = mExtRadius * sin ( phi );
245
246 vector point ( mPointEndCenterLine + XLocal * VecRadialDirection + YLocal *
VecRadialDirection2 ) ;
247 points [ NumCells + i] = point ;
248 phi += (1.0/ double ( NumCells )) *2* mPI ;
249 }
250
251
252 sprintf ( fileName ," actuatorDisk . vtk ");
253 file = fopen ( fileName ,"w");
254
255 fprintf ( file ,"# vtk DataFile Version 3.0\ n");
256 fprintf ( file ," Analytical surface of actuator disk . \n");
257 fprintf ( file ," ASCII \n");
258
259 fprintf ( file ," DATASET UNSTRUCTURED_GRID \n");
260 fprintf ( file ," POINTS %i float \n", NumPoints );
261
262 for ( int i = 0; i < points . size () ; i ++) {
263 fprintf ( file ,"%e %e %e\n", points [i ]. x () , points [i ]. y () , points [i ]. z () ) ;
264 }
265
266 fprintf ( file ," CELLS %i %i\n", NumCells , NumInts );
267
268 for ( unsigned int i = 0; i < NumCells -1; i ++) {
269 fprintf ( file ,"%i %i %i %i %i \n" ,4,i , i+ NumCells ,i+ NumCells +1 , i +1) ;
270 }
271 fprintf ( file ,"%i %i %i %i %i \n" ,4, NumCells -1 ,2* NumCells -1 , NumCells ,0) ;
272
273 fprintf ( file ," CELL_TYPES %i\n", NumCells );
274
275 for ( unsigned int i = 0; i < NumCells ; i ++) {
276 fprintf ( file ,"%i\n" ,9) ;
277 }
278
279 fclose ( file ) ;
280 }
281
282
283 bool actuatorDiskExplicitForce :: PointIsInDisk ( const vector & iPointStartCenterLine , const vector
& iPointEndCenterLine , const vector & iPoint , scalar & oDist2 , vector & oLineTangent , vector
& oCircumferentialDirection ) {
284 //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
285 // Check if a given point is located in the actuator disk region .
286 //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
287
288 vector VecLineTangent ( iPointEndCenterLine - iPointStartCenterLine );
289 scalar LineTangentLength = sqrt ( VecLineTangent .x () * VecLineTangent .x () + VecLineTangent .
y () * VecLineTangent .y () + VecLineTangent .z () * VecLineTangent . z () );
290
291 if( LineTangentLength != 0.0) {
292 VecLineTangent /= LineTangentLength ;
293 }
294 else {
295 Info << " Warning : The centerline tangent has zero length .\n";
296 return false ;
297 }
298
299 oLineTangent = VecLineTangent ;
300
301 vector VecStartLineToPoint ( iPoint - iPointStartCenterLine ) ;
302 scalar PointProjOnLine = VecStartLineToPoint & VecLineTangent ;
303
304 // Check if the point is inside the actuator disk in the axial direction
305 if (!( PointProjOnLine >= 0.0 && PointProjOnLine <= LineTangentLength )) {
306 return false ;
307 }
308
309 vector VecLineToPoint ( VecStartLineToPoint - ( VecLineTangent * PointProjOnLine ) );
310 scalar RadialDist2 = VecLineToPoint .x () * VecLineToPoint .x () + VecLineToPoint .y () *
VecLineToPoint .y () + VecLineToPoint .z () * VecLineToPoint . z () ;
311 oDist2 = RadialDist2 ;
312
313 oCircumferentialDirection = VecLineTangent ^ VecLineToPoint ;
314 oCircumferentialDirection /= mag ( oCircumferentialDirection );
315
316 // Check if the point is inside the actuator disk in the radial direction
317 return ( RadialDist2 <= mExtRadius * mExtRadius && RadialDist2 >= mIntRadius * mIntRadius );
318
319 }
320
321 bool actuatorDiskExplicitForce :: PointIsInHub ( const vector & iPointStartCenterLine , const vector
& iPointEndCenterLine , const vector & iPoint ) {
322 //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
323 // Check if a given point is located within the outer surface of the actuator disk region
and so close to the centerline
324 // that the radial distance is smaller than the interior radius of the actuator disk .
325 // This function is currently not used .
326 //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
327
328 vector VecLineTangent ( iPointEndCenterLine - iPointStartCenterLine );
329 scalar LineTangentLength = sqrt ( VecLineTangent .x () * VecLineTangent .x () + VecLineTangent .
y () * VecLineTangent .y () + VecLineTangent .z () * VecLineTangent . z () );
330
331 if( LineTangentLength != 0.0) {
332 VecLineTangent /= LineTangentLength ;
333 }
334 else {
335 Info << " Warning : The centerline tangent has zero length .\n";
336 return false ;
337 }
338
339 vector VecStartLineToPoint ( iPoint - iPointStartCenterLine ) ;
340 scalar PointProjOnLine = VecStartLineToPoint & VecLineTangent ;
341
342 // Check if the point is inside the actuator disk in the axial direction
343 if (!( PointProjOnLine >= 0.0 && PointProjOnLine <= LineTangentLength )) {
344 return false ;
345 }
346
347 vector VecLineToPoint ( VecStartLineToPoint - ( VecLineTangent * PointProjOnLine ) );
348 scalar RadialDist2 = VecLineToPoint .x () * VecLineToPoint .x () + VecLineToPoint .y () *
VecLineToPoint .y () + VecLineToPoint .z () * VecLineToPoint . z () ;
349
350 // Check if the point is inside the actuator disk in the radial direction
351 return ( RadialDist2 < mIntRadius * mIntRadius );
352
353 }
354
355
356 scalar actuatorDiskExplicitForce :: CalcAxialForce ( const scalar & iRadialDist , const scalar & iRho )
{
357 //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
358 // Compute the force component in the axial direction . The force is computed from a simple
equation resulting in a force
359 // that varies with the radial distance .
360 // If you have a better model of a rotor , comment the four lines below and add your own
calculation of the axial force .
361 // Do not forget to also change the calculation of the tangential force ( CalcCircForce ())
below .
362 //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
363 scalar axialForce = 0.0;
364 scalar radiusScaled = ( iRadialDist / mExtRadius - mIntRadius / mExtRadius ) /(1.0 -
mIntRadius / mExtRadius );
365 scalar Ax = (105.0/8.0) * mThrust /( CalcDiskThickness () * mPI *(3.0* mIntRadius +4.0* mExtRadius
) *( mExtRadius - mIntRadius ) );
366 axialForce = Ax * radiusScaled * sqrt (1.0 - radiusScaled );
367 //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
368
369 if( debug >= 2) {
370 Info << " Axial force : " << axialForce << "\n";
371 }
372
373 return axialForce ;
374 }
375
376 scalar actuatorDiskExplicitForce :: CalcCircForce ( const scalar & iRadialDist , const scalar & iRho )
{
377 //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
378 // Compute the force component in the tangential direction . The force is computed from a
simple equation resulting in a force
379 // that varies with the radial distance .
380 // Change the four lines below if you have a better model .
381 //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
382 scalar tangentialForce = 0.0;
383 scalar radiusScaled = ( iRadialDist / mExtRadius - mIntRadius / mExtRadius ) /(1.0 -
mIntRadius / mExtRadius );
384 scalar At = (105.0/8.0) * mTorque /( CalcDiskThickness () * mPI * mExtRadius *( mExtRadius -
mIntRadius ) *(3.0* mExtRadius +4.0* mIntRadius ));
385 tangentialForce = ( At * radiusScaled * sqrt (1.0 - radiusScaled ) /( radiusScaled *(1.0 -
mIntRadius / mExtRadius ) + mIntRadius / mExtRadius ) );
386 //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
387
388 if( debug >= 2) {
389 Info << " Tangential force : " << tangentialForce << "\n";
390 }
391
392 return tangentialForce ;
393 }
394
395
396 } // end namespace Foam
