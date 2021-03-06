#include "sample_points_to_grid.h"

#include "common/checks.h"

#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "Eigen/Dense"

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

vtkStandardNewMacro(sample_points_to_grid);

sample_points_to_grid::sample_points_to_grid()
{
    this->SetNumberOfInputPorts(2);
    this->SetNumberOfOutputPorts(1);
}

sample_points_to_grid::~sample_points_to_grid() {}

int sample_points_to_grid::FillInputPortInformation(int port, vtkInformation* info)
{
    if (!this->Superclass::FillInputPortInformation(port, info))
    {
        return 0;
    }

    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
        return 1;
    }
    if (port == 1)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkRectilinearGrid");
        return 1;
    }

    return 0;
}

int sample_points_to_grid::RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*)
{
    return 1;
}

int sample_points_to_grid::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get access to information and data
    auto in_points = input_vector[0]->GetInformationObject(0);
    auto points = vtkPointSet::SafeDownCast(in_points->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(points, "Input on port 0 is not a point set.");

    auto in_grid = input_vector[1]->GetInformationObject(0);
    auto grid = vtkRectilinearGrid::SafeDownCast(in_grid->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(grid, "Input on port 1 is not a rectilinear grid.");

    // Get output
    auto out_info = output_vector->GetInformationObject(0);
    auto output = vtkRectilinearGrid::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(output, "Output is not a rectilinear grid.");

    output->CopyStructure(grid);

    // Check parameters
    __check_param_min_ret(this->NumberOfBins[0], 1, "Number of bins must be at least 1.", 0);
    __check_param_min_ret(this->NumberOfBins[1], 1, "Number of bins must be at least 1.", 0);
    __check_param_min_ret(this->NumberOfBins[2], 1, "Number of bins must be at least 1.", 0);

    // Create coarse regular grid and use as for binning input points
    const auto x_num_bins = static_cast<long long>(this->NumberOfBins[0]) + 2;
    const auto y_num_bins = static_cast<long long>(this->NumberOfBins[1]) + 2;
    const auto z_num_bins = static_cast<long long>(this->NumberOfBins[2]) + 2;

    std::array<double, 6> bounds{};
    grid->GetBounds(bounds.data());

    const auto x_width = bounds[1] - bounds[0];
    const auto y_width = bounds[3] - bounds[2];
    const auto z_width = bounds[5] - bounds[4];

    const auto x_bin_width = x_width / (x_num_bins - 2);
    const auto y_bin_width = y_width / (y_num_bins - 2);
    const auto z_bin_width = z_width / (z_num_bins - 2);

    const auto x_origin = bounds[0] - x_bin_width;
    const auto y_origin = bounds[2] - y_bin_width;
    const auto z_origin = bounds[4] - z_bin_width;

    const auto x_upper_bound = bounds[1] + x_bin_width;
    const auto y_upper_bound = bounds[3] + y_bin_width;
    const auto z_upper_bound = bounds[5] + z_bin_width;

    const Eigen::Vector3d origin(x_origin, y_origin, z_origin);
    const Eigen::Vector3d bin_size(x_bin_width, y_bin_width, z_bin_width);

    // Create bins
    const auto grid_size = x_num_bins * y_num_bins * z_num_bins;

    std::vector<std::vector<vtkIdType>> bins(grid_size + 1);

    auto get_bin_index = [&](const Eigen::Vector3d& coords) -> std::tuple<vtkIdType, int, int, int>
    {
        // If points are out of bounds, return index of the dummy cell
        if (coords[0] < x_origin || coords[1] < y_origin || coords[2] < z_origin ||
            coords[0] > x_upper_bound || coords[1] > y_upper_bound || coords[2] > z_upper_bound)
        {
            return std::make_tuple(grid_size, -1, -1, -1);
        }

        // Calculate index from coordinates
        const Eigen::Vector3d cell_coords = (coords - origin).cwiseQuotient(bin_size);

        return std::make_tuple(static_cast<vtkIdType>(std::floor(cell_coords[0]) + x_num_bins * (std::floor(cell_coords[1]) + y_num_bins * std::floor(cell_coords[2]))),
            static_cast<int>(std::floor(cell_coords[0])), static_cast<int>(std::floor(cell_coords[1])), static_cast<int>(std::floor(cell_coords[2])));
    };

    // Create point data for all point/cell data of the input
    std::map<int, std::pair<vtkDataSetAttributes*, int>> array_map;

    for (int i = 0; i < points->GetPointData()->GetNumberOfArrays(); ++i)
    {
        const auto input_array = points->GetPointData()->GetArray(i);

        vtkSmartPointer<vtkDataArray> output_array;
        output_array.TakeReference(vtkDataArray::CreateDataArray(input_array->GetDataType()));
        output_array->SetName(input_array->GetName());
        output_array->SetNumberOfComponents(input_array->GetNumberOfComponents());
        output_array->SetNumberOfTuples(grid->GetNumberOfPoints());

        array_map[output->GetPointData()->AddArray(output_array)] = std::make_pair(points->GetPointData(), i);
    }

    for (int i = 0; i < points->GetCellData()->GetNumberOfArrays(); ++i)
    {
        const auto input_array = points->GetCellData()->GetArray(i);

        vtkSmartPointer<vtkDataArray> output_array;
        output_array.TakeReference(vtkDataArray::CreateDataArray(input_array->GetDataType()));
        output_array->SetName(input_array->GetName());
        output_array->SetNumberOfComponents(input_array->GetNumberOfComponents());
        output_array->SetNumberOfTuples(grid->GetNumberOfPoints());

        array_map[output->GetPointData()->AddArray(output_array)] = std::make_pair(points->GetCellData(), i);
    }

    // Create array storing the distance to the nearest neighbor as approximation for error
    auto error_array = vtkSmartPointer<vtkDoubleArray>::New();
    error_array->SetName("Error distance");
    error_array->SetNumberOfComponents(1);
    error_array->SetNumberOfTuples(grid->GetNumberOfPoints());

    // Iterate over points and assign it to a bin
    /**
        OpenMP information

        read global:    points, get_bin_index
        write global:   bins
    **/
    #pragma omp parallel for schedule(dynamic) default(none) shared(points, get_bin_index, bins)
    for (vtkIdType p = 0; p < points->GetNumberOfPoints(); ++p)
    {
        Eigen::Vector3d point;
        points->GetPoint(p, point.data());

        bins[std::get<0>(get_bin_index(point))].push_back(p);
    }

    // For each node of the grid, find closest point and use its value
    std::array<int, 6> extent{};
    grid->GetExtent(extent.data());

    auto x_coords = grid->GetXCoordinates();
    auto y_coords = grid->GetYCoordinates();
    auto z_coords = grid->GetZCoordinates();

    const auto num_x_nodes = extent[1] - extent[0] + 1;
    const auto num_y_nodes = extent[3] - extent[2] + 1;
    const auto num_z_nodes = extent[5] - extent[4] + 1;

    auto get_grid_index = [&](const vtkIdType i, const vtkIdType j, const vtkIdType k) -> vtkIdType
    {
        return static_cast<vtkIdType>(i + num_x_nodes * (j + num_y_nodes * k));
    };

    std::size_t num_bad_nodes = 0;

    /**
        OpenMP information

        read global:    array_map, bins, points, x_coords, y_coords, z_coords
        write global:   error_array, output
    **/
    #pragma omp parallel for schedule(dynamic) default(none) shared(error_array, output, \
        array_map, bins, points, x_coords, y_coords, z_coords)
    for (int k = extent[4]; k <= extent[5]; ++k)
    {
        const auto z = z_coords->GetComponent(k, 0);

        for (int j = extent[2]; j <= extent[3]; ++j)
        {
            const auto y = y_coords->GetComponent(j, 0);

            for (int i = extent[0]; i <= extent[1]; ++i)
            {
                const auto x = x_coords->GetComponent(i, 0);
                const Eigen::Vector3d coords(x, y, z);

                const auto bin_index = get_bin_index(coords);
                const auto grid_index = get_grid_index(i, j, k);

                if (std::get<0>(bin_index) != grid_size)
                {
                    double closest_point_distance = std::numeric_limits<double>::max();
                    vtkIdType closest_point = -1;

                    // Add neighboring bins to find closest point
                    std::vector<vtkIdType> cells{ std::get<0>(bin_index) };

                    for (int r = -1; r <= 1; ++r)
                    {
                        const auto z_offset = r * x_num_bins * y_num_bins;

                        if (std::get<3>(bin_index) + r >= 0 && static_cast<long long>(std::get<3>(bin_index)) + r < z_num_bins)
                        {
                            for (int q = -1; q <= 1; ++q)
                            {
                                const auto y_offset = q * x_num_bins;

                                if (std::get<2>(bin_index) + q >= 0 && static_cast<long long>(std::get<2>(bin_index)) + q < y_num_bins)
                                {
                                    for (int p = -1; p <= 1; ++p)
                                    {
                                        const auto x_offset = static_cast<long long>(p);

                                        if (std::get<1>(bin_index) + p >= 0 && static_cast<long long>(std::get<1>(bin_index)) + p < x_num_bins
                                            && (r != 0 || q != 0 || p != 0))
                                        {
                                            cells.push_back(std::get<0>(bin_index) + x_offset + y_offset + z_offset);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // Find closest point
                    Eigen::Vector3d point;

                    for (const auto& cell : cells)
                    {
                        for (const auto& p : bins[cell])
                        {
                            points->GetPoint(p, point.data());

                            if ((point - coords).norm() < closest_point_distance)
                            {
                                closest_point_distance = (point - coords).norm();
                                closest_point = p;
                            }
                        }
                    }

                    if (closest_point == -1)
                    {
                        #pragma omp critical(bad_nodes)
                        {
                            ++num_bad_nodes;
                        }

                        error_array->SetValue(grid_index, std::numeric_limits<double>::quiet_NaN());

                        continue;
                    }

                    // Copy properties of nearest neighbor
                    for (int a = 0; a < output->GetPointData()->GetNumberOfArrays(); ++a)
                    {
                        for (int c = 0; c < output->GetPointData()->GetArray(a)->GetNumberOfComponents(); ++c)
                        {
                            output->GetPointData()->GetArray(a)->SetComponent(grid_index, c,
                                array_map.at(a).first->GetArray(array_map.at(a).second)->GetComponent(closest_point, c));
                        }
                    }

                    error_array->SetValue(grid_index, closest_point_distance);
                }
            }
        }
    }

    if (num_bad_nodes > 0)
    {
        std::cerr << "Could not find closest points for " << num_bad_nodes << " nodes." << std::endl;
    }

    // Pass field data and add error array
    output->GetFieldData()->ShallowCopy(points->GetFieldData());
    output->GetPointData()->AddArray(error_array);

    return 1;
}
