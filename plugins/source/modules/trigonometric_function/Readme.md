# Trigonometric function

Polyline defined by a trigonometric function.

The sine function, e.g., is of the form ![Equation](https://render.githubusercontent.com/render/math?math={\alpha}*sin({\beta}*t%2B\gamma)%2B\delta).

## Parameters

The following parameters are available in the properties panel in ParaView:

| Parameter                 | Description                                                                           | Default value |
|---------------------------|---------------------------------------------------------------------------------------|---------------|
| Number of output points   | Number of points defining the output polyline.                                        | 100           |
| Length                    | Length of the center line.                                                            | 10            |
| Function                  | Enum offering different trigonometric functions.                                      | sine          |
| Alpha                     | Outer scalar factor, see equation above.                                              | 1             |
| Beta                      | Inner scalar factor, see equation above.                                              | 1             |
| Gamma                     | Inner offset, see equation above.                                                     | 0             |
| Delta                     | Outer offset, see equation above.                                                     | 0             |
