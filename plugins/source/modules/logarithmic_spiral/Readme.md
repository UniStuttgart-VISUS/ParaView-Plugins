# Trigonometric function

Polyline defined by the logarithmic spiral function.

The logarithmic spiral function is

(1) ![Equation](https://render.githubusercontent.com/render/math?math=r(\phi)=ae^{k\phi}),

with size factor ![Equation](https://render.githubusercontent.com/render/math?math=a) and spiral slope ![Equation](https://render.githubusercontent.com/render/math?math=k). The calculation of the cartesian coordinates is given by

(2) ![Equation](https://render.githubusercontent.com/render/math?math=x(\phi)=r(\phi)\\,\\,\text{cos}\\,\phi), and ![Equation](https://render.githubusercontent.com/render/math?math=y(\phi)=r(\phi)\\,\\,\text{sin}\\,\phi).

## Parameters

The following parameters are available in the properties panel in ParaView:

| Parameter                 | Description                                                                           | Default value |
|---------------------------|---------------------------------------------------------------------------------------|---------------|
| Number of output points   | Number of points defining the output polyline.                                        | 100           |
| Length                    | Length of the spiral.                                                                 | 10            |
| Size factor               | Size factor multiplied with the spiral function.                                      | 1             |
| Slope                     | Slope of the spiral.                                                                  | 1             |
