# B-Spline

Create a B-spline from the input de Boor points.

## Input

The following inputs can be connected to the filter:

| Input                     | Description                                                                               | Type          | Remark        |
|---------------------------|-------------------------------------------------------------------------------------------|---------------|---------------|
| Points                    | De Boor points defining the B-spline.                                                     | Poly data     |               |
| Point                     | Point for which the closest corresponding point on the B-spline should be calculated.     | Poly data     | optional      |

## Parameters

The following parameters are available in the properties panel in ParaView:

| Parameter                 | Description                                                                                                   | Default value         |
|---------------------------|---------------------------------------------------------------------------------------------------------------|-----------------------|
| Number of output points   | Number of output points defining the B-spline polyline.                                                       | 100                   |
| Degree                    | Degree of the B-spline.                                                                                       | 2                     |
| Hit endpoints             | Use multiplicity for reaching the first and last de Boor point.                                               | yes                   |
