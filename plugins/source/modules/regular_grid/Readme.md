# Regular Grid

Creates a regular grid with given extent, origin and cell sizes. Its main purpose is to use for resampling other data onto this grid.

## Parameters

The following parameters can be set by the user:

| Parameter     | Description                                       | Type          | Accepted values       | Default value     |
| ------------- | ------------------------------------------------- | ------------- | -----------------     | ----------------- |
| Output type   | Output grid type.                                 | Enumeration   | Regular, Rectilinear  | Regular           |
| Anchor        | Defines where the origin of the grid is located.  | Enumeration   | Lowest, Center        | Lowest            |
| Dimension     | Number of cells in the respective direction.      | Integer       | > 0                   | \[16, 16, 16\]    |
| Origin        | Origin of the grid.                               | Float         |                       | \[0, 0, 0\]       |
| Cell size     | Cell size in the respective direction.            | Float         | > 0                   | \[1, 1, 1\]       |

The anchor set to 'lowest' coresponds to the origin of the grid being at the smallest coordinate (lower left back). Contrary 'center' indicates that the grid's center coresponds to the origin.

## Output

The output is a uniform grid, stored as either a regular (vtkImageData) or rectilinear grid (vtkRectilinearGrid).
