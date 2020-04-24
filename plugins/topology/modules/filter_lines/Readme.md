# Filter lines

Filter lines using different criteria.

## Input

The following inputs can be connected to the filter:

| Input                     | Description                                                               | Type          | Remark        |
|---------------------------|---------------------------------------------------------------------------|---------------|---------------|
| Lines                     | Poly lines which should be filtered.                                      | Poly data     |               |

## Parameters

The following parameters are available in the properties panel in ParaView:

| Parameter                 | Description                                                                           | Default value |
|---------------------------|---------------------------------------------------------------------------------------|---------------|
| AngleFilter               | Filter by angle between the line segment and a user-defined vector array.             | off           |
| MaxAngle                  | Maximum angle for which the segment is considered valid.                              | 45            |
| AngleVector               | Name of the array storing the vectors used for the angle calculation.                 |               |
| SizeFilter                | Filter polylines by number of segments.                                               | on            |
| MinSize                   | Minumum number of segments for which a polyline is considered valid                   | 5             |
