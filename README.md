# ParaView Plugins

ParaView plugins for small, recurrent tasks.

## List of plugins

The sources and filters are grouped into separate plugins, based on their functionality:

| Plugin                            | Description                                                               |
|-----------------------------------|---------------------------------------------------------------------------|
| [data](#data-plugin)              | Data loading, saving and transformation.                                  |
| [geometry](#geometry-plugin)      | Geometrical computations and modifications.                               |
| [source](#source-plugin)          | Simple data sources for geometry and test dataset.                        |
| [topology](#topology-plugin)      | Computation and extraction of topological features.                       |

### Sources and filters

The following plugins provide sources and filters.

#### Data plugin

This plugin provides filters for data loading and saving, as well for transformation between data types.

| Filter                                                                                                        | Description                                                                           |
|---------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------|
| [Grid to polydata](plugins/data/modules/grid_to_polydata/Readme.md)                                           | Convert an unstructured grid to polydata.                                             |
| [Image to rectilinear grid](plugins/data/modules/image_to_rectilinear_grid/Readme.md)                         | Convert an image to a rectilinear grid.                                               |
| [Rectilinear grid to structured grid](plugins/data/modules/rectilinear_grid_to_structured_grid/Readme.md)     | Convert an rectilinear grid to a structured grid.                                     |
| [Resample rotating grid](plugins/data/modules/resample_rotating_grid/Readme.md)                               | Given a rotation around the z-axis, the grid is resampled as if it did not rotate.    |
| [Sample points to grid](plugins/data/modules/sample_points_to_grid/Readme.md)                                 | Using nearest neighbor interpolation, sample points to grid nodes.                    |
| [Separate block](plugins/data/modules/separate_block/Readme.md)                                               | Extract a separate block in its native data structure.                                |
| [Structured grid to unstructured grid](plugins/data/modules/structured_grid_to_unstructured_grid/Readme.md)   | Convert a structured grid to an unstructured grid.                                    |

#### Geometry plugin

This plugin provides filters for the computation and modification of geometrical structures, such as lines.

| Filter                                                                            | Description                                                                           |
|-----------------------------------------------------------------------------------|---------------------------------------------------------------------------------------|
| [B-Spline](plugins/geometry/modules/b_spline/Readme.md)                           | Create a B-spline from the input de Boor points.                                      |
| [Connect lines](plugins/geometry/modules/connect_lines/Readme.md)                 | Connect single lines into polylines, where points are shared.                         |
| [Connect points](plugins/geometry/modules/connect_points/Readme.md)               | Connect points into a polyline, in the order they are stored.                         |
| [Open closed lines](plugins/geometry/modules/open_closed_lines/Readme.md)         | Open closed lines by removing the last segment or by duplicating the shared point.    |
| [Point cells](plugins/geometry/modules/point_cells/Readme.md)                     | Points for which a cell should be created.                                            |
| [Remove long segments](plugins/geometry/modules/remove_long_segments/Readme.md)   | Remove segments from a polyline, which are by far larger than the median segment.     |
| [Smooth lines](plugins/geometry/modules/smooth_lines/Readme.md)                   | Smooth lines using Gaussian or Taubin smoothing.                                      |
| [Sort line points](plugins/geometry/modules/sort_line_points/Readme.md)           | Sort stored points, thus they appear in the order of occurance in the lines.          |
| [Truncate lines](plugins/geometry/modules/truncate_lines/Readme.md)               | Truncate lines by defining an offset and a number of points for the new lines.        |

#### Source plugin

This plugin provides sources for geometry and test datasets.

| Source                                                                            | Description                                                                                   |
|-----------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| [Helix](plugins/source/modules/helix/Readme.md)                                   | Polyline in form of a simple helix with user-defined length, radius and number of windings.   |
| [Logarithmic spiral](plugins/source/modules/logarithmic_spiral/Readme.md)         | Polyline defined by the logarithmic spiral function.                                          |
| [Regular grid](plugins/source/modules/regular_grid/Readme.md)                     | Empty regular grid, defined by number of cells, origin and cell sizes.                        |
| [Trigonometric function](plugins/source/modules/trigonometric_function/Readme.md) | Polyline defined by a trigonometric function.                                                 |

#### Topology plugin

This plugin provides filters for the computation and extraction of topological features, mainly for vector fields of flow simulations.

| Filter                                                                                        | Description                                                                                                           |
|-----------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------|
| [Attachment/Separation lines](plugins/topology/modules/attachment_separation_lines/Readme.md) | Provide a scalar field on a two-dimensional surface, which zero-contour indicates attachment and separation lines.    |
| [Filter lines](plugins/topology/modules/filter_lines/Readme.md)                               | Filter lines using differenc criteria, such as size and angle.                                                        |

## Usage

The project uses CMake as build system. To configure the project run CMake. The following configuration options are available:

| Option                                | Description                                           | Default value     |
|---------------------------------------|-------------------------------------------------------|-------------------|
| PARAVIEW_PLUGIN_ENABLE_VISUSdata      | Enable the [data](#data-plugin) plugin.               | on                |
| PARAVIEW_PLUGIN_ENABLE_VISUSgeometry  | Enable the [geometry](#geometry-plugin) plugin.       | on                |
| PARAVIEW_PLUGIN_ENABLE_VISUSsource    | Enable the [source](#source-plugin) plugin.           | on                |
| PARAVIEW_PLUGIN_ENABLE_VISUStopology  | Enable the [topology](#topology-plugin) plugin.       | on                |

# License

This project is published under the MIT license. In the following you can find a list of contributors.

## List of contributors

- Alexander Straub, University of Stuttgart, Germany  
  (alexander.straub@visus.uni-stuttgart.de)

## MIT license

Copyright (c) 2020 University of Stuttgart Visualization Research Center

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
