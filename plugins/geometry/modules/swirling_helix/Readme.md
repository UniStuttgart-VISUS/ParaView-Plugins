# Swirling helix

Create a helix that swirls around an input B-spline of degree &gt;2.

## Input

The following inputs can be connected to the filter:

| Input                     | Description                                                               | Type          | Remark        |
|---------------------------|---------------------------------------------------------------------------|---------------|---------------|
| Lines                     | Lines around the created helix should swirl.                              | Poly data     |               |

## Parameters

The following parameters are available in the properties panel in ParaView:

| Parameter                 | Description                                                                                                   | Default value         |
|---------------------------|---------------------------------------------------------------------------------------------------------------|-----------------------|
| Tangent                   | Tangents of the input B-spline.                                                                               |                       |
| Normal                    | Normals of the input B-spline.                                                                                |                       |
| Radius                    | Radius, i.e., distance of the swirling helix to the input line.                                               | 1                     |
| Windings                  | Number of windings of the swirling helix around the input line.                                               | 2                     |
