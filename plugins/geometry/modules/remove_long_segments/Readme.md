# Remove long segments

Remove segments that are longer by a factor than the median segment of a polyline. Removing a segment results in two individual, unconnected lines.

## Input

The following inputs can be connected to the filter:

| Input                     | Description                                                               | Type          | Remark        |
|---------------------------|---------------------------------------------------------------------------|---------------|---------------|
| Lines                     | Lines from which segments are removed.                                    | Poly data     |               |

## Parameters

The following parameters are available in the properties panel in ParaView:

| Parameter         | Description                                                                                                   | Default value         |
|-------------------|---------------------------------------------------------------------------------------------------------------|-----------------------|
| Length factor     | Factor determining the removal of a segment, if it is larger than the factor times the median segment length. | 2.0                   |
