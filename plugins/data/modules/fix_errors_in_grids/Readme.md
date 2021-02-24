# Fix errors in grids

Using an error estimation, fix values at nodes exceeding an error threshold by linear interpolation.

## Input

The following inputs can be connected to the filter:

| Input                     | Description                                                                               | Type                  | Remark        |
|---------------------------|-------------------------------------------------------------------------------------------|-----------------------|---------------|
| Input                     | Grid with error estimation, whose nodes should be fixed.                                  | Rectilinear grid      |               |

## Parameters

The following parameters are available in the properties panel in ParaView:

| Parameter                 | Description                                                                                                   | Default value         |
|---------------------------|---------------------------------------------------------------------------------------------------------------|-----------------------|
| Errors                    | Array containing the error estimation.                                                                        |                       |
