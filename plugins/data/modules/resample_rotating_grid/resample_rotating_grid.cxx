#include "resample_rotating_grid.h"

#include "common/checks.h"
#include "common/pointers.h"

#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkCommand.h"
#include "vtkDataArray.h"
#include "vtkDataArraySelection.h"
#include "vtkFieldData.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTable.h"

#include "Eigen/Dense"

#include <array>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

vtkStandardNewMacro(resample_rotating_grid);

resample_rotating_grid::resample_rotating_grid()
{
    this->SetNumberOfInputPorts(2);
    this->SetNumberOfOutputPorts(1);

    this->RotationAngle = nullptr;
    this->RotationalVelocity = nullptr;

    this->scalar_array_selection = vtkSmartPointer<vtkDataArraySelection>::New();
    this->scalar_array_selection->AddObserver(vtkCommand::ModifiedEvent, this, &vtkObject::Modified);

    this->vector_array_selection = vtkSmartPointer<vtkDataArraySelection>::New();
    this->vector_array_selection->AddObserver(vtkCommand::ModifiedEvent, this, &vtkObject::Modified);

    this->velocity_array_selection = vtkSmartPointer<vtkDataArraySelection>::New();
    this->velocity_array_selection->AddObserver(vtkCommand::ModifiedEvent, this, &vtkObject::Modified);

    this->pass_point_array_selection = vtkSmartPointer<vtkDataArraySelection>::New();
    this->pass_point_array_selection->AddObserver(vtkCommand::ModifiedEvent, this, &vtkObject::Modified);

    this->pass_cell_array_selection = vtkSmartPointer<vtkDataArraySelection>::New();
    this->pass_cell_array_selection->AddObserver(vtkCommand::ModifiedEvent, this, &vtkObject::Modified);
}

resample_rotating_grid::~resample_rotating_grid()
{
    __delete_ptr(this->RotationAngle);
    __delete_ptr(this->RotationalVelocity);
}

vtkDataArraySelection* resample_rotating_grid::GetScalarArraySelection()
{
    return this->scalar_array_selection;
}

vtkDataArraySelection* resample_rotating_grid::GetVectorArraySelection()
{
    return this->vector_array_selection;
}

vtkDataArraySelection* resample_rotating_grid::GetVelocityArraySelection()
{
    return this->velocity_array_selection;
}

vtkDataArraySelection* resample_rotating_grid::GetPassPointArraySelection()
{
    return this->pass_point_array_selection;
}

vtkDataArraySelection* resample_rotating_grid::GetPassCellArraySelection()
{
    return this->pass_cell_array_selection;
}

int resample_rotating_grid::FillInputPortInformation(int port, vtkInformation* info)
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
    if (port == 1)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
        return 1;
    }

    return 0;
}

int resample_rotating_grid::RequestUpdateExtent(vtkInformation* vtkNotUsed(request), vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    return 1;
}

int resample_rotating_grid::RequestData(vtkInformation* vtkNotUsed(request), vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get access to information and data
    auto in_info = input_vector[0]->GetInformationObject(0);
    auto input = vtkRectilinearGrid::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(input, "Input on port 0 is not a rectilinear grid.");

    auto in_table = input_vector[1]->GetInformationObject(0);
    auto table = vtkTable::SafeDownCast(in_table->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(table, "Input on port 1 is not a table.");

    // Get output
    auto out_info = output_vector->GetInformationObject(0);
    auto output = vtkRectilinearGrid::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(output, "Output is not a rectilinear grid.");

    output->CopyStructure(input);

    // Create array indicating validity of the node
    auto validity = vtkSmartPointer<vtkIdTypeArray>::New();
    validity->SetName("Validity");
    validity->SetNumberOfComponents(1);
    validity->SetNumberOfTuples(input->GetNumberOfPoints());

    // Get user-selected data arrays
    std::vector<std::pair<vtkSmartPointer<vtkDataArray>, vtkDataArray*>> scalars, vectors, velocities;
    std::vector<vtkSmartPointer<vtkDataArray>> point_pass_through, cell_pass_through, field_pass_through;

    for (int index = 0; index < input->GetPointData()->GetNumberOfArrays(); ++index)
    {
        auto in_array = input->GetPointData()->GetArray(index);

        if (in_array && in_array->GetName())
        {
            if (this->scalar_array_selection->ArrayIsEnabled(in_array->GetName()))
            {
                vtkSmartPointer<vtkDataArray> scalar_array;
                scalar_array.TakeReference(vtkDataArray::CreateDataArray(in_array->GetDataType()));
                scalar_array->DeepCopy(in_array);

                scalars.push_back(std::make_pair(scalar_array, in_array));
            }
            else if (this->vector_array_selection->ArrayIsEnabled(in_array->GetName()))
            {
                vtkSmartPointer<vtkDataArray> vector_array;
                vector_array.TakeReference(vtkDataArray::CreateDataArray(in_array->GetDataType()));
                vector_array->DeepCopy(in_array);

                vectors.push_back(std::make_pair(vector_array, in_array));
            }
            else if (this->velocity_array_selection->ArrayIsEnabled(in_array->GetName()))
            {
                vtkSmartPointer<vtkDataArray> velocity_array;
                velocity_array.TakeReference(vtkDataArray::CreateDataArray(in_array->GetDataType()));
                velocity_array->DeepCopy(in_array);

                velocities.push_back(std::make_pair(velocity_array, in_array));
            }
            else if (this->pass_point_array_selection->ArrayIsEnabled(in_array->GetName()))
            {
                vtkSmartPointer<vtkDataArray> pass_array;
                pass_array.TakeReference(vtkDataArray::CreateDataArray(in_array->GetDataType()));
                pass_array->DeepCopy(in_array);

                point_pass_through.push_back(pass_array);
            }
        }
    }

    for (int index = 0; index < input->GetCellData()->GetNumberOfArrays(); ++index)
    {
        auto in_array = input->GetCellData()->GetArray(index);

        if (in_array && in_array->GetName())
        {
            if (this->pass_cell_array_selection->ArrayIsEnabled(in_array->GetName()))
            {
                vtkSmartPointer<vtkDataArray> pass_array;
                pass_array.TakeReference(vtkDataArray::CreateDataArray(in_array->GetDataType()));
                pass_array->DeepCopy(in_array);

                cell_pass_through.push_back(pass_array);
            }
        }
    }

    for (int index = 0; index < input->GetFieldData()->GetNumberOfArrays(); ++index)
    {
        auto in_array = input->GetFieldData()->GetArray(index);

        if (in_array && in_array->GetName())
        {
            vtkSmartPointer<vtkDataArray> pass_array;
            pass_array.TakeReference(vtkDataArray::CreateDataArray(in_array->GetDataType()));
            pass_array->DeepCopy(in_array);

            field_pass_through.push_back(pass_array);
        }
    }

    // Get rotation information from table
    __check_not_null_ret(this->RotationAngle, "Column not set for rotation angle.");
    __check_not_null_ret(this->RotationalVelocity, "Column not set for angular frequency.");

    auto rotation_column = table->GetColumnByName(this->RotationAngle);
    auto frequency_column = table->GetColumnByName(this->RotationalVelocity);

    __check_not_null_ret(rotation_column, "Column not found for rotation angle.");
    __check_not_null_ret(frequency_column, "Column not found for angular frequency.");

    __check_not_empty_ret(rotation_column, "Column for rotation angle does not have an entry.", 0);
    __check_not_empty_ret(frequency_column, "Column for angular frequency does not have an entry.", 0);

    const auto angle = rotation_column->GetVariantValue(0).ToDouble();
    const auto omega = frequency_column->GetVariantValue(0).ToDouble();

    // Create rotation matrix
    Eigen::Matrix3d rotation, inverse_rotation;

    rotation <<
        std::cos(angle), -std::sin(angle), 0.0,
        std::sin(angle), std::cos(angle), 0.0,
        0.0, 0.0, 1.0;

    inverse_rotation = rotation.inverse();

    // Get nodes of output grid as points
    std::array<int, 6> extent{};
    output->GetExtent(extent.data());

    auto x_coords = output->GetXCoordinates();
    auto y_coords = output->GetYCoordinates();
    auto z_coords = output->GetZCoordinates();

    for (int k = extent[4]; k <= extent[5]; ++k)
    {
        for (int j = extent[2]; j <= extent[3]; ++j)
        {
            for (int i = extent[0]; i <= extent[1]; ++i)
            {
                const auto index = i + ((extent[1] - extent[0] + 1) * (j + (extent[3] - extent[2] + 1) * k));

                // Rotate coordinates
                const Eigen::Vector3d coords(
                    x_coords->GetComponent(i, 0),
                    y_coords->GetComponent(j, 0),
                    z_coords->GetComponent(k, 0)
                );

                const Eigen::Vector3d rotated_coords = rotation * coords;

                // Variables for interpolation
                int sub_id;
                std::array<double, 3> p_coords{};
                std::array<double, 8> weights{};

                const auto cell = input->FindAndGetCell(const_cast<Eigen::Vector3d&>(rotated_coords).data(),
                    nullptr, 0, 0.00001, sub_id, p_coords.data(), weights.data());

                if (cell != nullptr)
                {
                    validity->SetValue(index, 1);

                    // Resample scalar arrays
                    for (auto& scalar_array : scalars)
                    {
                        scalar_array.first->InterpolateTuple(index, cell->GetPointIds(), scalar_array.second, weights.data());
                    }

                    // Resample vector arrays
                    for (auto& vector_array : vectors)
                    {
                        __check_num_components_ret(vector_array.first, 3, "Vectors must be three-dimensional.");

                        vector_array.first->InterpolateTuple(index, cell->GetPointIds(), vector_array.second, weights.data());

                        // Apply Jacobian of the transformation
                        Eigen::Vector3d vec(
                            vector_array.first->GetComponent(index, 0),
                            vector_array.first->GetComponent(index, 1),
                            vector_array.first->GetComponent(index, 2));

                        vec = inverse_rotation * vec;

                        vector_array.first->SetComponent(index, 0, vec[0]);
                        vector_array.first->SetComponent(index, 1, vec[1]);
                        vector_array.first->SetComponent(index, 2, vec[2]);
                    }

                    // Resample velocities
                    for (auto& velocity_array : velocities)
                    {
                        __check_num_components_ret(velocity_array.first, 3, "Velocities must be three-dimensional.");

                        velocity_array.first->InterpolateTuple(index, cell->GetPointIds(), velocity_array.second, weights.data());

                        // Remove rotational velocity part and then apply Jacobian of the transformation
                        Eigen::Vector3d vec(
                            velocity_array.first->GetComponent(index, 0),
                            velocity_array.first->GetComponent(index, 1),
                            velocity_array.first->GetComponent(index, 2));

                        vec -= omega * Eigen::Vector3d(-rotated_coords[1], rotated_coords[0], 0.0);
                        vec = inverse_rotation * vec;

                        velocity_array.first->SetComponent(index, 0, vec[0]);
                        velocity_array.first->SetComponent(index, 1, vec[1]);
                        velocity_array.first->SetComponent(index, 2, vec[2]);
                    }
                }
                else
                {
                    validity->SetValue(index, 0);

                    // Set zero
                    for (auto& scalar_array : scalars)
                    {
                        scalar_array.first->SetComponent(index, 0, 0.0);
                    }

                    for (auto& vector_array : vectors)
                    {
                        vector_array.first->SetComponent(index, 0, 0.0);
                        vector_array.first->SetComponent(index, 1, 0.0);
                        vector_array.first->SetComponent(index, 2, 0.0);
                    }

                    for (auto& velocity_array : velocities)
                    {
                        velocity_array.first->SetComponent(index, 0, 0.0);
                        velocity_array.first->SetComponent(index, 1, 0.0);
                        velocity_array.first->SetComponent(index, 2, 0.0);
                    }
                }
            }
        }
    }

    // Add arrays to output
    for (auto& scalar_array : scalars)
    {
        output->GetPointData()->AddArray(scalar_array.first);
    }

    for (auto& vector_array : vectors)
    {
        output->GetPointData()->AddArray(vector_array.first);
    }

    for (auto& velocity_array : velocities)
    {
        output->GetPointData()->AddArray(velocity_array.first);
    }

    for (auto& pass_through_array : point_pass_through)
    {
        output->GetPointData()->AddArray(pass_through_array);
    }

    for (auto& pass_through_array : cell_pass_through)
    {
        output->GetCellData()->AddArray(pass_through_array);
    }

    for (auto& pass_through_array : field_pass_through)
    {
        output->GetFieldData()->AddArray(pass_through_array);
    }

    output->GetPointData()->AddArray(validity);

    return 1;
}
