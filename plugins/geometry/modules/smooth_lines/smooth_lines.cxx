#include "smooth_lines.h"

#include "vtkDataObject.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStructuredGrid.h"

#include "Eigen/Dense"

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

vtkStandardNewMacro(smooth_lines);

smooth_lines::smooth_lines()
{
    this->SetNumberOfInputPorts(2);
    this->SetNumberOfOutputPorts(1);
}

int smooth_lines::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
        return 1;
    }
    else if (port == 1)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkStructuredGrid");
        info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
        return 1;
    }

    return 0;
}

int smooth_lines::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    return 1;
}

int smooth_lines::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get input
    auto input = vtkPolyData::GetData(input_vector[0]);

    vtkSmartPointer<vtkStructuredGrid> perp_grid = nullptr;

    if (input_vector[1] != nullptr && vtkStructuredGrid::GetData(input_vector[1]) != nullptr)
    {
        auto grid = vtkStructuredGrid::GetData(input_vector[1]);

        auto vector_field = GetInputArrayToProcess(0, grid);

        if (vector_field != nullptr)
        {
            auto perp_vector_field = vtkSmartPointer<vtkDoubleArray>::New();
            perp_vector_field->DeepCopy(vector_field);

            // Rotate vectors 90 degrees
            Eigen::Vector3d temp_vec;
            double temp{};

            for (vtkIdType i = 0; i < perp_vector_field->GetNumberOfTuples(); ++i)
            {
                perp_vector_field->GetTuple(i, temp_vec.data());

                temp = temp_vec[0];
                temp_vec[0] = temp_vec[1];
                temp_vec[1] = -temp;

                perp_vector_field->SetTuple(i, temp_vec.data());
            }

            perp_grid = vtkSmartPointer<vtkStructuredGrid>::New();
            perp_grid->CopyStructure(grid);
            perp_grid->GetPointData()->AddArray(perp_vector_field);
        }
    }

    // Get output
    auto* out_info = output_vector->GetInformationObject(0);
    auto* output = vtkPolyData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    // Check parameters
    if (this->NumIterations < 0)
    {
        std::cerr << "Number of iterations must not be negative" << std::endl;
        return 0;
    }

    const auto method = static_cast<method_t>(this->Method);

    if (this->Lambda <= 0.0 || this->Lambda > 1.0)
    {
        std::cerr << "The smoothing factor must be in the range (0, 1]" << std::endl;
        return 0;
    }

    if (method == method_t::taubin && (this->Mu < -1.0 || this->Mu > 0.0))
    {
        std::cerr << "The inflation factor must be in the range [-1, 0]" << std::endl;
        return 0;
    }

    if (method == method_t::taubin && (std::abs(this->Mu) <= this->Lambda && this->Mu != 0.0))
    {
        std::cerr << "The absolute value of the inflation factor must be larger than the smoothing factor, or equal 0" << std::endl;
        return 0;
    }

    if (method == method_t::perp_gaussian && perp_grid == nullptr)
    {
        std::cerr << "Input vector field must be available for perpendicular Gaussian smoothing methods" << std::endl;
        return 0;
    }

    // Copy input to output and work on that copy
    output->DeepCopy(input);

    // Line-wise
    output->GetLines()->InitTraversal();

    auto point_indices = vtkSmartPointer<vtkIdList>::New();

    while (output->GetLines()->GetNextCell(point_indices))
    {
        std::vector<Eigen::Vector3d> points(point_indices->GetNumberOfIds()), temp_points(point_indices->GetNumberOfIds());

        for (vtkIdType index = 0; index < point_indices->GetNumberOfIds(); ++index)
        {
            output->GetPoints()->GetPoint(point_indices->GetId(index), points[index].data());
        }

        // Apply smoothing...
        for (std::size_t i = 0; i < this->NumIterations; ++i)
        {
            // ... Gaussian smoothing
            if (method == method_t::gaussian || method == method_t::taubin)
            {
                for (std::size_t index = 0; index < points.size(); ++index)
                {
                    temp_points[index] = gaussian_smoothing(points, index, this->Lambda);
                }
            }
            // ... perpendicular Gaussian
            else if (method == method_t::perp_gaussian)
            {
                for (std::size_t index = 0; index < points.size(); ++index)
                {
                    temp_points[index] = gaussian_smoothing(points, index, this->Lambda, perp_grid);
                }
            }

            // ... second step for Taubin smoothing
            if (method == method_t::taubin)
            {
                for (std::size_t index = 0; index < points.size(); ++index)
                {
                    points[index] = gaussian_smoothing(temp_points, index, this->Mu);
                }
            }
            else
            {
                std::swap(points, temp_points);
            }
        }

        // Write smoothed points to output
        for (vtkIdType index = 0; index < point_indices->GetNumberOfIds(); ++index)
        {
            output->GetPoints()->SetPoint(point_indices->GetId(index), points[index].data());
        }
    }

    return 1;
}

Eigen::Vector3d smooth_lines::gaussian_smoothing(const std::vector<Eigen::Vector3d>& points,
    const std::size_t index, const double weight) const
{
    Eigen::Vector3d direction;

    if (index == 0)
    {
        direction = weight * (points[index + 1] - points[index]);
    }
    else if (index == points.size() - 1)
    {
        direction = weight * (points[index - 1] - points[index]);
    }
    else
    {
        direction = weight * (
            0.5 * (points[index - 1] - points[index]) +
            0.5 * (points[index + 1] - points[index]));
    }

    return points[index] + direction;
}

Eigen::Vector3d smooth_lines::gaussian_smoothing(const std::vector<Eigen::Vector3d>& points,
    const std::size_t index, const double weight, const vtkStructuredGrid* vector_field) const
{
    if (index == 0 || index == points.size() - 1)
    {
        return points[index];
    }

    Eigen::Vector3d direction = gaussian_smoothing(points, index, weight) - points[index];

    // Interpolate vector field at the original point position
    int sub_id{};
    std::array<double, 3> p_coords{};
    std::vector<double> weights(10); // support cells with up to 10 points for now

    auto* cell = const_cast<vtkStructuredGrid*>(vector_field)->FindAndGetCell(const_cast<double*>(points[index].data()),
        nullptr, 0, 1e-10, sub_id, p_coords.data(), weights.data());

    assert(cell->GetNumberOfPoints() <= 10);

    Eigen::Vector3d vector, vector_interpolated;
    vector_interpolated.setZero();

    for (vtkIdType i = 0; i < cell->GetNumberOfPoints(); ++i)
    {
        const_cast<vtkStructuredGrid*>(vector_field)->GetPointData()->GetArray(0)->GetTuple(cell->GetPointId(i), vector.data());

        vector_interpolated += weights[i] * vector;
    }

    vector_interpolated.normalize();

    // Project direction onto vector field to restrict movement
    direction = direction.dot(vector_interpolated) * vector_interpolated;

    return points[index] + direction;
}
