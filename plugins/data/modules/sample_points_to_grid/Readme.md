# Resample rotating grid

Given the angle of rotation around the z-axis, the grid is resampled as if it did not rotate.

## Input

The following inputs can be connected to the filter:

| Input                     | Description                                                                               | Type                  | Remark        |
|---------------------------|-------------------------------------------------------------------------------------------|-----------------------|---------------|
| Input                     | Input grid, which should be resampled.                                                    | Rectilinear grid      |               |
| Rotation                  | Table containing the information about the rotation.                                      | Table                 |               |

## Parameters

The following parameters are available in the properties panel in ParaView:

| Parameter                 | Description                                                                                                   | Default value         |
|---------------------------|---------------------------------------------------------------------------------------------------------------|-----------------------|
| Rotation (angle)          | Column containing the angle of rotation.                                                                      |                       |
| Scalar arrays             | Selection of scalar arrays to resample.                                                                       |                       |
| Vector arrays             | Selection of vector arrays to resample.                                                                       |                       |
| Pass point arrays         | Selection of point arrays to simply pass through without modification.                                        |                       |
| Pass cell arrays          | Selection of cell arrays to simply pass through without modification.                                         |                       |
