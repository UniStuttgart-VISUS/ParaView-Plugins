# ParaView Plugins

ParaView plugins for small, recurrent tasks.

## List of plugins

The sources and filters are grouped into separate plugins, based on their functionality:

| Plugin                            | Description                                                               |
|-----------------------------------|---------------------------------------------------------------------------|
| [source](#source-plugin)          | Simple data sources for geometry and test dataset.                        |
| [topology](#topology-plugin)      | Computation and extraction of topological features.                       |

### Sources and filters

The following plugins provide sources and filters.

#### Source plugin

This plugin provides sources for geometry and test datasets.

| Source                                            | Description                                                                                  |
|---------------------------------------------------|----------------------------------------------------------------------------------------------|
| [Helix](plugins/source/modules/helix/Readme.md)   | Polyline in form of a simple helix with user-defined length, radius and number of windings.  |

#### Topology plugin

This plugin provides filters for the computation and extraction of topological features, mainly for vector fields of flow simulations.

| Source                                                                                        | Description                                                                                                           |
|-----------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------|
| [Attachment/Separation lines](plugins/topology/modules/attachment_separation_lines/Readme.md) | Provide a scalar field on a two-dimensional surface, which zero-contour indicates attachment and separation lines.    |

## Usage

The project uses CMake as build system. To configure the project run CMake. The following configuration options are available:

| Option                                | Description                                           | Default value     |
|---------------------------------------|-------------------------------------------------------|-------------------|
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
