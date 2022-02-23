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
#include "Eigen/Sparse"

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
    const auto variant = static_cast<variant_t>(this->Variant);

    if (this->ImplicitLambda <= 0.0)
    {
        std::cerr << "The implicit smoothing factor must be positive" << std::endl;
        return 0;
    }

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
    vtkSmartPointer<vtkPolyData> output_steps;

    {
        output_steps = vtkSmartPointer<vtkPolyData>::New();
        output_steps->DeepCopy(input);

        auto original_positions = vtkSmartPointer<vtkDoubleArray>::New();
        original_positions->SetName("Original Positions");
        original_positions->SetNumberOfComponents(3);
        original_positions->SetNumberOfTuples(input->GetNumberOfPoints());

        std::array<double, 3> point{};

        for (vtkIdType p = 0; p < input->GetNumberOfPoints(); ++p)
        {
            input->GetPoint(p, point.data());
            original_positions->SetTuple(p, point.data());
        }

        auto displacements = vtkSmartPointer<vtkDoubleArray>::New();
        displacements->SetName("Displacements");
        displacements->SetNumberOfComponents(3);
        displacements->SetNumberOfTuples(input->GetNumberOfPoints());

        output_steps->GetPointData()->AddArray(original_positions);
        output_steps->GetPointData()->AddArray(displacements);
    }

    // Line-wise
    output_steps->GetLines()->InitTraversal();

    auto point_indices = vtkSmartPointer<vtkIdList>::New();

    while (output_steps->GetLines()->GetNextCell(point_indices))
    {
        std::vector<Eigen::Vector3d> points(point_indices->GetNumberOfIds()), temp_points(point_indices->GetNumberOfIds());

        for (vtkIdType index = 0; index < point_indices->GetNumberOfIds(); ++index)
        {
            output_steps->GetPoints()->GetPoint(point_indices->GetId(index), points[index].data());
        }

        variant_t state = variant;
        auto max_distance = (points[0] - points[points.size() - 1]).norm();

        if (method == method_t::implicit_gaussian)
        {
            // Create matrix representing the points of the polyline
            const auto n = points.size();

            Eigen::MatrixXd vertices(n, 3);

            for (std::size_t i = 0; i < n; ++i)
            {
                vertices.row(i) = points[i].transpose();
            }

            // Create weight matrix for fixed end points
            Eigen::SparseMatrix<double> L_fixed(n, n);

            for (std::size_t j = 1; j < n - 1; ++j)
            {
                const auto weight_left = 1.0 / (points[j] - points[j - 1]).norm();
                const auto weight_right = 1.0 / (points[j + 1] - points[j]).norm();
                const auto weight_sum = weight_left + weight_right;

                L_fixed.insert(j, j - 1) = weight_left / weight_sum;
                L_fixed.insert(j, j) = -1.0;
                L_fixed.insert(j, j + 1) = weight_right / weight_sum;
            }

            // Create weight matrix for moving end points
            auto L_moving = L_fixed;

            L_moving.insert(0, 0) = -1.0;
            L_moving.insert(0, 1) = 1.0;
            L_moving.insert(n - 1, n - 2) = 1.0;
            L_moving.insert(n - 1, n - 1) = -1.0;

            // Create matrices
            Eigen::SparseMatrix<double> I(n, n);
            I.setIdentity();

            const Eigen::SparseMatrix<double> A_fixed = (I - this->ImplicitLambda * L_fixed);
            const Eigen::SparseMatrix<double> A_moving = (I - this->ImplicitLambda * L_moving);

            // Apply smoothing
            for (std::size_t i = 0; i < this->NumIterations; ++i)
            {
                switch (state)
                {
                case variant_t::fixed_endpoints:
                    gaussian_smoothing(vertices, A_fixed);

                    break;
                case variant_t::growing:
                case variant_t::normal:
                    gaussian_smoothing(vertices, A_moving);
                }

                if (state == variant_t::growing)
                {
                    const auto distance = (vertices.row(0) - vertices.row(points.size() - 1)).norm();

                    if (distance < 0.9 * max_distance)
                    {
                        state = variant_t::fixed_endpoints;
                    }

                    max_distance = std::max(max_distance, distance);
                }
            }

            // Write smoothed points to output
            auto original_positions = input->GetPointData()->GetArray("Original Positions");
            auto displacements = output_steps->GetPointData()->GetArray("Displacements");

            for (vtkIdType index = 0; index < point_indices->GetNumberOfIds(); ++index)
            {
                Eigen::Vector3d original_position;
                original_positions->GetTuple(point_indices->GetId(index), original_position.data());

                const Eigen::Vector3d point = vertices.row(index).transpose();
                const Eigen::Vector3d displacement = point - original_position;

                output_steps->GetPoints()->SetPoint(point_indices->GetId(index), point.data());

                displacements->SetTuple(point_indices->GetId(index), displacement.data());
            }
        }
        else
        {
            // Apply smoothing...
            for (std::size_t i = 0; i < this->NumIterations; ++i)
            {
                // ... Gaussian smoothing
                if (method == method_t::gaussian || method == method_t::taubin)
                {
                    for (std::size_t index = 0; index < points.size(); ++index)
                    {
                        temp_points[index] = gaussian_smoothing(points, index, this->Lambda, state);
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
                        points[index] = gaussian_smoothing(temp_points, index, this->Mu, state);
                    }
                }
                else
                {
                    std::swap(points, temp_points);
                }

                if (state == variant_t::growing)
                {
                    const auto distance = (points[0] - points[points.size() - 1]).norm();

                    if (distance < 0.9 * max_distance)
                    {
                        state = variant_t::fixed_endpoints;
                    }

                    max_distance = std::max(max_distance, distance);
                }
            }

            // Write smoothed points to output
            auto original_positions = input->GetPointData()->GetArray("Original Positions");
            auto displacements = output_steps->GetPointData()->GetArray("Displacements");

            for (vtkIdType index = 0; index < point_indices->GetNumberOfIds(); ++index)
            {
                Eigen::Vector3d original_position;
                original_positions->GetTuple(point_indices->GetId(index), original_position.data());

                const Eigen::Vector3d point = points[index];
                const Eigen::Vector3d displacement = point - original_position;

                output_steps->GetPoints()->SetPoint(point_indices->GetId(index), point.data());

                displacements->SetTuple(point_indices->GetId(index), displacement.data());
            }
        }
    }

    return 1;
}

void smooth_lines::gaussian_smoothing(Eigen::MatrixXd& vertices,
    const Eigen::SparseMatrix<double>& matrix) const
{
    const Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver(matrix);

    const Eigen::VectorXd x1 = solver.solve(vertices.col(0));
    const Eigen::VectorXd x2 = solver.solve(vertices.col(1));
    const Eigen::VectorXd x3 = solver.solve(vertices.col(2));

    vertices << x1, x2, x3;
}

Eigen::Vector3d smooth_lines::gaussian_smoothing(const std::vector<Eigen::Vector3d>& points,
    const std::size_t index, const double weight, const variant_t variant) const
{
    Eigen::Vector3d direction;

    if (index == 0)
    {
        if (variant == variant_t::fixed_endpoints)
        {
            direction.setZero();
        }
        else
        {
            direction = weight * (points[index + 1] - points[index]);
        }
    }
    else if (index == points.size() - 1)
    {
        if (variant == variant_t::fixed_endpoints)
        {
            direction.setZero();
        }
        else
        {
            direction = weight * (points[index - 1] - points[index]);
        }
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

    Eigen::Vector3d direction = gaussian_smoothing(points, index, weight, variant_t::fixed_endpoints) - points[index];

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
