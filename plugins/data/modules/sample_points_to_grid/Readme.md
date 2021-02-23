# Sample points to grid

Using nearest neighbor interpolation, sample points to grid nodes.

## Input

The following inputs can be connected to the filter:

| Input                     | Description                                                                               | Type                  | Remark        |
|---------------------------|-------------------------------------------------------------------------------------------|-----------------------|---------------|
| Input                     | Point input whose data is resampled on the grid.                                          | Point set             |               |
| Grid                      | Grid at whose nodes the points are resampled.                                             | Rectilinear grid      |               |

## Parameters

The following parameters are available in the properties panel in ParaView:

| Parameter                 | Description                                                                                                   | Default value         |
|---------------------------|---------------------------------------------------------------------------------------------------------------|-----------------------|
| Number of bins            | Number of bins per spatial dimension for nearest neighbor search, reducing complexity.                        | 10                    |
