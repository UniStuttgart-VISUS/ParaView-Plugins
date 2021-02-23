#include "fix_errors_in_grids.h"

#include "common/checks.h"

#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "Eigen/Dense"

#include <algorithm>
#include <array>
#include <limits>
#include <set>

vtkStandardNewMacro(fix_errors_in_grids);

fix_errors_in_grids::fix_errors_in_grids()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

fix_errors_in_grids::~fix_errors_in_grids() {}

int fix_errors_in_grids::FillInputPortInformation(int port, vtkInformation* info)
{
    if (!this->Superclass::FillInputPortInformation(port, info))
    {
        return 0;
    }

    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkRectilinearGrid");
        return 1;
    }

    return 0;
}

int fix_errors_in_grids::RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*)
{
    return 1;
}

int fix_errors_in_grids::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get access to information and data
    auto in_grid = input_vector[0]->GetInformationObject(0);
    auto grid = vtkRectilinearGrid::SafeDownCast(in_grid->Get(vtkDataObject::DATA_OBJECT()));

    // Get output
    auto out_info = output_vector->GetInformationObject(0);
    auto output = vtkRectilinearGrid::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    output->DeepCopy(grid);

    // Get array containing error estimation
    auto errors = GetInputArrayToProcess(0, grid);

    __check_not_null_ret(errors, "Array containing error estimation not selected.");
    __check_num_components_ret(errors, 1, "Array containing error estimation must have one component.");

    // Iterate over nodes and calculate minimum and median error values
    double min_error, max_error, median_error;
    min_error = std::numeric_limits<double>::max();
    max_error = 0.0;

    std::multiset<double> error_values;

    std::size_t num_values = 0;

    for (vtkIdType i = 0; i < errors->GetNumberOfValues(); ++i)
    {
        const auto error = errors->GetComponent(i, 0);

        if (error >= 0.0)
        {
            min_error = std::min(min_error, error);
            max_error = std::max(max_error, error);
            error_values.insert(error);
        }
    }

    const auto median_it = std::next(error_values.begin(), std::distance(error_values.begin(), error_values.end()) / 2);

    if (median_it != error_values.end())
    {
        median_error = *median_it;
    }
    else
    {
        median_error = *error_values.begin();
    }

    // Set threshold for "good" values as triple the difference between minimum and median
    const auto error_threshold = min_error + 3 * (median_error - min_error);

    // Correct values
    std::array<int, 6> extent;
    grid->GetExtent(extent.data());

    auto x_coords = grid->GetXCoordinates();
    auto y_coords = grid->GetYCoordinates();
    auto z_coords = grid->GetZCoordinates();

    const auto num_x_nodes = static_cast<long long>(extent[1]) - extent[0] + 1;
    const auto num_y_nodes = static_cast<long long>(extent[3]) - extent[2] + 1;
    const auto num_z_nodes = static_cast<long long>(extent[5]) - extent[4] + 1;

    std::size_t num_corrected_nodes = 0;

    for (vtkIdType k = extent[4]; k <= extent[5]; ++k)
    {
        const auto z = z_coords->GetComponent(k, 0);

        for (vtkIdType j = extent[2]; j <= extent[3]; ++j)
        {
            const auto y = y_coords->GetComponent(j, 0);

            for (vtkIdType i = extent[0]; i <= extent[1]; ++i)
            {
                const auto x = x_coords->GetComponent(i, 0);
                const Eigen::Vector3d coords(x, y, z);

                const auto grid_index = i + num_x_nodes * (j + num_y_nodes * k);

                // Check error and compare to threshold
                const auto error = errors->GetComponent(grid_index, 0);

                if (error > error_threshold && !(i == extent[0] || j == extent[2] || k == extent[4] || i == extent[1] || j == extent[3] || k == extent[5]))
                {
                    // Get neighbors and compute barycenter of triangle in each adjacent cell
                    const auto index_left = grid_index - 1;
                    const auto index_right = grid_index + 1;
                    const auto index_bottom = grid_index - num_x_nodes;
                    const auto index_top = grid_index + num_x_nodes;
                    const auto index_back = grid_index - num_x_nodes * num_y_nodes;
                    const auto index_front = grid_index + num_x_nodes * num_y_nodes;

                    const auto weight_left = (x_coords->GetComponent(i + 1, 0) - x) / (x_coords->GetComponent(i + 1, 0) - x_coords->GetComponent(i - 1, 0));
                    const auto weight_right = (x - x_coords->GetComponent(i - 1, 0)) / (x_coords->GetComponent(i + 1, 0) - x_coords->GetComponent(i - 1, 0));
                    const auto weight_bottom = (y_coords->GetComponent(j + 1, 0) - y) / (y_coords->GetComponent(j + 1, 0) - y_coords->GetComponent(j - 1, 0));
                    const auto weight_top = (y - y_coords->GetComponent(j - 1, 0)) / (y_coords->GetComponent(j + 1, 0) - y_coords->GetComponent(j - 1, 0));
                    const auto weight_back = (z_coords->GetComponent(k + 1, 0) - z) / (z_coords->GetComponent(k + 1, 0) - z_coords->GetComponent(k - 1, 0));
                    const auto weight_front = (z - z_coords->GetComponent(k - 1, 0)) / (z_coords->GetComponent(k + 1, 0) - z_coords->GetComponent(k - 1, 0));

                    for (int pd = 0; pd < grid->GetPointData()->GetNumberOfArrays(); ++pd)
                    {
                        const auto in_array = grid->GetPointData()->GetArray(pd);
                        const auto out_array = output->GetPointData()->GetArray(pd);
                        const auto num_components = in_array->GetNumberOfComponents();

                        Eigen::Matrix<double, Eigen::Dynamic, 1> left(num_components);
                        Eigen::Matrix<double, Eigen::Dynamic, 1> right(num_components);
                        Eigen::Matrix<double, Eigen::Dynamic, 1> bottom(num_components);
                        Eigen::Matrix<double, Eigen::Dynamic, 1> top(num_components);
                        Eigen::Matrix<double, Eigen::Dynamic, 1> back(num_components);
                        Eigen::Matrix<double, Eigen::Dynamic, 1> front(num_components);

                        in_array->GetTuple(index_left, left.data());
                        in_array->GetTuple(index_right, right.data());
                        in_array->GetTuple(index_bottom, bottom.data());
                        in_array->GetTuple(index_top, top.data());
                        in_array->GetTuple(index_back, back.data());
                        in_array->GetTuple(index_front, front.data());

                        const Eigen::Matrix<double, Eigen::Dynamic, 1> interpolated = (1.0 / 3.0) * (
                            weight_left * left + weight_right * right +
                            weight_bottom * bottom + weight_top * top +
                            weight_back * back + weight_front * front);

                        out_array->SetTuple(grid_index, interpolated.data());
                    }

                    // Set dummy error value to indicate that the node has been repaired
                    output->GetPointData()->GetArray(errors->GetName())->SetComponent(grid_index, 0, -1.0);

                    ++num_corrected_nodes;
                }
            }
        }
    }

    return 1;
}
