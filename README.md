# viscoelasticFlow
Porting the viscoelastic solver from foam-extend-4.1 to OpenFOAM-v1912.

Out-dated code and deprecated features are replaced with the latest standard.

## Instructions:
1. Download the source code, (not neccesary to save in OpenFOAM folder).
2. Compile the viscoelastic library. Run `wmake` in `viscoelasticModels` folder.
3. Compile the viscoelastic solver. Run `wmake` in `viscoelasticFluidFoam` folder.
4. Now you can test tutorial cases in `tutorials`

## Acknowledgement
I would like to thank Jovani L.Favero who wrote the original solver.
