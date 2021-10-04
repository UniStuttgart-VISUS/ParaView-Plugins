#include "image_to_rectilinear_grid.h"

#include "common/checks.h"

#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <array>

vtkStandardNewMacro(image_to_rectilinear_grid);

image_to_rectilinear_grid::image_to_rectilinear_grid()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

image_to_rectilinear_grid::~image_to_rectilinear_grid() {}

int image_to_rectilinear_grid::FillInputPortInformation(int port, vtkInformation* info)
{
    if (!this->Superclass::FillInputPortInformation(port, info))
    {
        return 0;
    }

    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
        return 1;
    }

    return 0;
}

int image_to_rectilinear_grid::RequestUpdateExtent(vtkInformation* vtkNotUsed(request), vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    return 1;
}

int image_to_rectilinear_grid::RequestData(vtkInformation* vtkNotUsed(request), vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get access to information and data
    auto in_info = input_vector[0]->GetInformationObject(0);
    auto input = vtkImageData::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(input, "Input is not an image.");

    // Get output
    auto out_info = output_vector->GetInformationObject(0);
    auto output = vtkRectilinearGrid::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(output, "Output is not a rectilinear grid.");

    // Create rectilinear grid from image
    std::array<int, 6> extent{};
    std::array<double, 3> origin{}, spacing{};

    input->GetExtent(extent.data());
    input->GetOrigin(origin.data());
    input->GetSpacing(spacing.data());

    auto coordinates_x = vtkSmartPointer<vtkDoubleArray>::New();
    auto coordinates_y = vtkSmartPointer<vtkDoubleArray>::New();
    auto coordinates_z = vtkSmartPointer<vtkDoubleArray>::New();

    coordinates_x->SetNumberOfComponents(1);
    coordinates_y->SetNumberOfComponents(1);
    coordinates_z->SetNumberOfComponents(1);

    coordinates_x->SetNumberOfTuples(static_cast<vtkIdType>(extent[1]) - extent[0] + 1);
    coordinates_y->SetNumberOfTuples(static_cast<vtkIdType>(extent[3]) - extent[2] + 1);
    coordinates_z->SetNumberOfTuples(static_cast<vtkIdType>(extent[5]) - extent[4] + 1);

    for (vtkIdType i = 0; i < coordinates_x->GetNumberOfTuples(); ++i)
    {
        coordinates_x->SetValue(i, origin[0] + i * spacing[0]);
    }

    for (vtkIdType i = 0; i < coordinates_y->GetNumberOfTuples(); ++i)
    {
        coordinates_y->SetValue(i, origin[1] + i * spacing[1]);
    }

    for (vtkIdType i = 0; i < coordinates_z->GetNumberOfTuples(); ++i)
    {
        coordinates_z->SetValue(i, origin[2] + i * spacing[2]);
    }

    output->SetExtent(extent.data());
    output->SetXCoordinates(coordinates_x);
    output->SetYCoordinates(coordinates_y);
    output->SetZCoordinates(coordinates_z);

    // Copy data
    output->GetPointData()->ShallowCopy(input->GetPointData());
    output->GetCellData()->ShallowCopy(input->GetCellData());
    output->GetFieldData()->ShallowCopy(input->GetFieldData());

    return 1;
}
