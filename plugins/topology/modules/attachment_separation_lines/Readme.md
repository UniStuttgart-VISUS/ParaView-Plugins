# Attachment and separation lines

This filter allows to extract attachment and separation lines on two-dimensional surfaces. This is done by projecting the velocity and the Jacobian matrix onto the surface and computing the determinant of the matrix

(1) ![Equation](https://render.githubusercontent.com/render/math?math=\left(\mathbf{u}_{2D}\,\,\,\,\,J_{2D}\mathbf{u}_{2D}\right)).

The resulting value is stored in a scalar field, where the zero-contour coincides with the attachment and separation lines.

## Input

The following inputs can be connected to the filter:

| Input                     | Description                                                               | Type          | Remark        |
|---------------------------|---------------------------------------------------------------------------|---------------|---------------|
| Surface mesh              | Surface mesh providing the normals, velocities and Jacobian matrices.     | Poly data     |               |

## Parameters

The following parameters are available in the properties panel in ParaView:

| Parameter                 | Description                                                                           | Default value |
|---------------------------|---------------------------------------------------------------------------------------|---------------|
| Normal                    | Name of the array storing the surface normals.                                        |               |
| Velocity                  | Name of the array storing the velocities.                                             |               |
| Jacobian                  | Name of the array storing the Jacobian matrices.                                      |               |
