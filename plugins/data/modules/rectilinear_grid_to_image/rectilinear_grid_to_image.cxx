#include "rectilinear_grid_to_image.h"

#include "common/checks.h"

#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkFieldData.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <array>
#include <cmath>
#include <iostream>

vtkStandardNewMacro(rectilinear_grid_to_image);

rectilinear_grid_to_image::rectilinear_grid_to_image()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

rectilinear_grid_to_image::~rectilinear_grid_to_image() {}

int rectilinear_grid_to_image::FillInputPortInformation(int port, vtkInformation* info)
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

int rectilinear_grid_to_image::RequestUpdateExtent(vtkInformation* vtkNotUsed(request), vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    return 1;
}

int rectilinear_grid_to_image::RequestData(vtkInformation* vtkNotUsed(request), vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get access to information and data
    auto in_info = input_vector[0]->GetInformationObject(0);
    auto input = vtkRectilinearGrid::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(input, "Input is not a rectilinear grid.");

    // Get output
    auto out_info = output_vector->GetInformationObject(0);
    auto output = vtkImageData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(output, "Output is not an image.");

    // Create image from rectilinear grid
    std::array<int, 6> extent{};

    input->GetExtent(extent.data());

    auto coordinates_x = input->GetXCoordinates();
    auto coordinates_y = input->GetYCoordinates();
    auto coordinates_z = input->GetZCoordinates();

    const std::array<double, 3> origin{
        coordinates_x->GetComponent(0, 0),
        coordinates_y->GetComponent(0, 0),
        coordinates_z->GetComponent(0, 0)
    };

    const std::array<double, 3> spacing{
        (coordinates_x->GetNumberOfTuples() > 1) ? (coordinates_x->GetComponent(1, 0) - coordinates_x->GetComponent(0, 0)) : 0.0,
        (coordinates_y->GetNumberOfTuples() > 1) ? (coordinates_y->GetComponent(1, 0) - coordinates_y->GetComponent(0, 0)) : 0.0,
        (coordinates_z->GetNumberOfTuples() > 1) ? (coordinates_z->GetComponent(1, 0) - coordinates_z->GetComponent(0, 0)) : 0.0
    };

    output->SetExtent(extent.data());
    output->SetOrigin(origin.data());
    output->SetSpacing(spacing.data());

    // Check if input is an image
    double previous = coordinates_x->GetComponent(0, 0);

    for (vtkIdType i = 1; i < coordinates_x->GetNumberOfTuples(); ++i)
    {
        const auto current = coordinates_x->GetComponent(i, 0);
        const auto difference = std::abs(current - previous);

        if (std::abs(spacing[0] - difference) > 1.0e-6)
        {
            std::cerr << "Input cannot be converted to image: X coordinates are not equally spaced." << std::endl;

            return 0;
        }

        previous = current;
    }

    previous = coordinates_y->GetComponent(0, 0);

    for (vtkIdType i = 1; i < coordinates_y->GetNumberOfTuples(); ++i)
    {
        const auto current = coordinates_y->GetComponent(i, 0);
        const auto difference = std::abs(current - previous);

        if (std::abs(spacing[0] - difference) > 1.0e-6)
        {
            std::cerr << "Input cannot be converted to image: Y coordinates are not equally spaced." << std::endl;

            return 0;
        }

        previous = current;
    }

    previous = coordinates_z->GetComponent(0, 0);

    for (vtkIdType i = 1; i < coordinates_z->GetNumberOfTuples(); ++i)
    {
        const auto current = coordinates_z->GetComponent(i, 0);
        const auto difference = std::abs(current - previous);

        if (std::abs(spacing[0] - difference) > 1.0e-6)
        {
            std::cerr << "Input cannot be converted to image: Z coordinates are not equally spaced." << std::endl;

            return 0;
        }

        previous = current;
    }

    // Copy data
    output->GetPointData()->ShallowCopy(input->GetPointData());
    output->GetCellData()->ShallowCopy(input->GetCellData());
    output->GetFieldData()->ShallowCopy(input->GetFieldData());

    return 1;
}
