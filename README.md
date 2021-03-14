# viscoelasticFlow
Porting the viscoelastic solver from foam-extend-4.1 to OpenFOAM-v2012.

Out-dated code and deprecated features are replaced with new methods.

Eigen v3.3.9 is included in the source code.

## Instructions:
1. Download the source code, (not neccesary to save in OpenFOAM folder).
2. Compile the viscoelastic library. Run `wmake` in `viscoelasticModels` folder.
3. Compile the viscoelastic solver. Run `wmake` in `viscoelasticFluidFoam` folder.
4. Now you can test tutorial cases in `tutorials`

## Acknowledgement
I would like to thank Jovani L.Favero who wrote the original solver.
