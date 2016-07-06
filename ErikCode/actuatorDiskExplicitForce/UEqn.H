    // Momentum predictor

//    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
/*
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
*/
    );

    UEqn().relax();

//Added by me
//Clear  old  values  in  ioVolumeForce
for(label k = 0; k < VolumeForce.size(); k++) {
VolumeForce[k] = vector::zero;
}

// Calculate  volume  force  from  actuator  disk
actuatorDisk.CalcActuatorDiskVolForce(mesh, VolumeForce);

fvOptions.constrain(UEqn());

solve(UEqn() == -fvc::grad(p) + VolumeForce);

fvOptions.correct(U);
//Finished by me


/* Original Code
    fvOptions.constrain(UEqn());

    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
*/
