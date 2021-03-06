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

#include "actuatorDiskExplicitForce.H"

#include "faceAreaPairGAMGAgglomeration.H"
#include "fvMesh.H"
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Write the outer surface of the actuator disk to a VTK file so that it can be visualized in Paraview .
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	FILE *file;
	char fileName[100];

	unsigned int NumCells = 20; // The cylindrical surface is visualized as 20 rectangular surfaces
	unsigned int NumPoints = 40; // 40 points are required ; 20 points at each end of the cylinder
	unsigned int NumInts = 5*NumCells; // Number of integers needed in the the VTK file ; each surface has 4 corner points , so we need 4 corner indices + the index of the surface = 5 indices per surface

	vectorField points(NumPoints, vector::zero);

	vector VecLineTangent(mPointEndCenterLine - mPointStartCenterLine);
	scalar LineTangentLength = sqrt(VecLineTangent.x()*VecLineTangent.x() + VecLineTangent.y()*VecLineTangent.y() + VecLineTangent.z()*VecLineTangent.z());

	if(LineTangentLength!= 0.0) {
		VecLineTangent /= LineTangentLength;
	}
	else {
		Info << "Warning: The centerline tangent has zero length.\n";
		return;	
	}

	// We need to find a vector in the radial direction . This can be any vector as long as it points in the radial direction .
	// First try with (1 0 0) and see if we can project it onto the normal plane of the actuator disk resulting in a vector in
	// the radial direction .
	vector VecRadialDirection (1.0,0.0,0.0);
	VecRadialDirection -= (VecRadialDirection & VecLineTangent)*VecLineTangent;

	if(mag(VecRadialDirection) < SMALL) {
		// If we enter this if statement , our guess (1 0 0) was parallel to the centerline of the actuator disk . Then
		// we try (0 1 0) instead . Since (1 0 0) was parallel to the centerline, (0 1 0) will for sure not be parallel to
		// the centerline .
		VecRadialDirection.x() = 0.0;
		VecRadialDirection.y() = 1.0;
		VecRadialDirection.z() = 0.0;

		VecRadialDirection -= (VecRadialDirection & VecLineTangent)*VecLineTangent;
	}

	if(mag(VecRadialDirection) > SMALL) {
		VecRadialDirection /= mag(VecRadialDirection);
	}
	else {
		Info << "Warning in actuatorDiskExplicitForce::WriteVTK(): mag(VecRadialDirection) close to zero.\n";
	}

	vector VecRadialDirection2 = VecLineTangent ^ VecRadialDirection;
	scalar XLocal = 0.0, YLocal = 0.0;

	// Compute points on first side of disk region
	double phi = 0.0;
	for(unsigned int i = 0; i < NumCells; i ++) {
		XLocal = mExtRadius*cos(phi);
		YLocal = mExtRadius*sin(phi);
	
		vector point(mPointStartCenterLine + XLocal*VecRadialDirection + YLocal*VecRadialDirection2);
		points[i] = point;
		phi += (1.0/ double(NumCells))*2*mPI ;
	}

	// Compute points on second side of disk region
	phi = 0.0;
	for(unsigned int i = 0; i < NumCells; i++) {
		XLocal = mExtRadius*cos(phi);
		YLocal = mExtRadius*sin(phi);
	
		vector point(mPointEndCenterLine + XLocal*VecRadialDirection + YLocal*VecRadialDirection2);
		points[NumCells + i] = point;
		phi += (1.0/ double(NumCells))*2*mPI;
	}


	sprintf(fileName, "actuatorDisk.vtk");
	file = fopen(fileName,"w");

	fprintf(file, "# vtk DataFile Version 3.0\n");
	fprintf(file, "Analytical surface of actuator disk.\n");
	fprintf(file, "ASCII\n");

	fprintf(file, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(file, "POINTS %i float\n", NumPoints);

	for(int i = 0; i < points.size(); i++) {
	fprintf(file, "%e %e %e\n", points[i].x(), points[i].y(), points[i].z()) ;
	}

	fprintf(file, "CELLS %i %i\n", NumCells, NumInts);
	
	for(unsigned int i = 0; i < NumCells-1; i++) {
		fprintf(file, "%i %i %i %i %i\n",4,i,i+NumCells,i+NumCells+1,i+1);
	}
	fprintf(file,"%i %i %i %i %i\n",4,NumCells-1,2*NumCells-1,NumCells,0);

	fprintf(file,"CELL_TYPES %i\n", NumCells);

	for(unsigned int i = 0; i < NumCells; i++) {
		fprintf(file, "%i\n",9);
	}
	
	fclose(file);
}


bool actuatorDiskExplicitForce::PointIsInDisk(const vector &iPointStartCenterLine, const vector &iPointEndCenterLine, const vector &iPoint, scalar &oDist2, vector &oLineTangent, vector &oCircumferentialDirection) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Check if a given point is located in the actuator disk region .
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	vector VecLineTangent(iPointEndCenterLine - iPointStartCenterLine);
	scalar LineTangentLength = sqrt(VecLineTangent.x()*VecLineTangent.x() + VecLineTangent.y()*VecLineTangent.y() + VecLineTangent.z()*VecLineTangent.z());

	if(LineTangentLength != 0.0) {
		VecLineTangent /= LineTangentLength;
	}
	else {
		Info << "Warning: The centerline tangent has zero length.\n";
		return false;
	}

	oLineTangent = VecLineTangent;

	vector VecStartLineToPoint(iPoint - iPointStartCenterLine);
	scalar PointProjOnLine = VecStartLineToPoint & VecLineTangent;

	// Check if the point is inside the actuator disk in the axial direction
	if(!(PointProjOnLine >= 0.0 && PointProjOnLine <= LineTangentLength)) {
		return false;
	}

	vector VecLineToPoint(VecStartLineToPoint - (VecLineTangent*PointProjOnLine));
	scalar RadialDist2 = VecLineToPoint.x()*VecLineToPoint.x() + VecLineToPoint.y()*VecLineToPoint.y() + VecLineToPoint.z()*VecLineToPoint.z();
	oDist2 = RadialDist2;

	oCircumferentialDirection = VecLineTangent ^ VecLineToPoint;
	oCircumferentialDirection /= mag (oCircumferentialDirection);
	
	// Check if the point is inside the actuator disk in the radial direction
	return(RadialDist2 <= mExtRadius*mExtRadius && RadialDist2 >= mIntRadius*mIntRadius);
		
}

bool actuatorDiskExplicitForce::PointIsInHub(const vector &iPointStartCenterLine, const vector &iPointEndCenterLine, const vector &iPoint) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Check if a given point is located within the outer surface of the actuator disk region and so close to the centerline
// that the radial distance is smaller than the interior radius of the actuator disk .
// This function is currently not used .
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	vector VecLineTangent(iPointEndCenterLine - iPointStartCenterLine);
	scalar LineTangentLength = sqrt(VecLineTangent.x()*VecLineTangent.x() + VecLineTangent.y()*VecLineTangent.y() + VecLineTangent.z()*VecLineTangent.z());
	
	if(LineTangentLength != 0.0) {
		VecLineTangent /= LineTangentLength;
	}
	else {
		Info << "Warning: The centerline tangent has zero length.\n";
	return false;
	}
	
	vector VecStartLineToPoint(iPoint - iPointStartCenterLine);
	scalar PointProjOnLine = VecStartLineToPoint & VecLineTangent;

	// Check if the point is inside the actuator disk in the axial direction
	if(!(PointProjOnLine >= 0.0 && PointProjOnLine <= LineTangentLength)) {
		return false;
	}

	vector VecLineToPoint(VecStartLineToPoint - (VecLineTangent*PointProjOnLine));
	scalar RadialDist2 = VecLineToPoint.x()*VecLineToPoint.x() + VecLineToPoint.y()*VecLineToPoint.y() + VecLineToPoint.z()*VecLineToPoint.z();

	// Check if the point is inside the actuator disk in the radial direction
	return(RadialDist2 < mIntRadius*mIntRadius);
	
}


scalar actuatorDiskExplicitForce::CalcAxialForce(const scalar &iRadialDist, const scalar &iRho) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute the force component in the axial direction . The force is computed from a simple equation resulting in a force
// that varies with the radial distance .
// If you have a better model of a rotor , comment the four lines below and add your own calculation of the axial force .
// Do not forget to also change the calculation of the tangential force ( CalcCircForce ()) below .
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	scalar axialForce = 0.0;
	scalar radiusScaled = (iRadialDist/mExtRadius - mIntRadius/mExtRadius)/(1.0 - mIntRadius/mExtRadius);
	scalar Ax = (105.0/8.0)*mThrust/( CalcDiskThickness()*mPI*(3.0*mIntRadius + 4.0*mExtRadius)*(mExtRadius - mIntRadius));
	axialForce = Ax*radiusScaled*sqrt(1.0 - radiusScaled);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(debug >= 2) {
		Info << "Axial force: " << axialForce << "\n";
	}

	return axialForce;
}

scalar actuatorDiskExplicitForce::CalcCircForce(const scalar &iRadialDist, const scalar &iRho) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute the force component in the tangential direction . The force is computed from a simple equation resulting in a force
// that varies with the radial distance .
// Change the four lines below if you have a better model .
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	scalar tangentialForce = 0.0;
	scalar radiusScaled = (iRadialDist/mExtRadius - mIntRadius/mExtRadius)/(1.0 - mIntRadius/mExtRadius);
	scalar At = (105.0/8.0)*mTorque/(CalcDiskThickness()*mPI*mExtRadius*(mExtRadius - mIntRadius)*(3.0*mExtRadius + 4.0*mIntRadius));
	tangentialForce = (At*radiusScaled*sqrt(1.0 - radiusScaled )/(radiusScaled*(1.0 - mIntRadius/mExtRadius) + mIntRadius/mExtRadius));
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	if(debug >= 2) {
		Info << "Tangential force: " << tangentialForce << "\n";
	}

	return tangentialForce;
}


} // end namespace Foam
