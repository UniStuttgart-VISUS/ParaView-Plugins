#include "rectilinear_grid_to_structured_grid.h"

#include "common/checks.h"

#include "vtkCellData.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStructuredGrid.h"

#include <array>

vtkStandardNewMacro(rectilinear_grid_to_structured_grid);

rectilinear_grid_to_structured_grid::rectilinear_grid_to_structured_grid()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

rectilinear_grid_to_structured_grid::~rectilinear_grid_to_structured_grid() {}

int rectilinear_grid_to_structured_grid::FillInputPortInformation(int port, vtkInformation* info)
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

int rectilinear_grid_to_structured_grid::RequestUpdateExtent(vtkInformation* vtkNotUsed(request), vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    return 1;
}

int rectilinear_grid_to_structured_grid::RequestData(vtkInformation* vtkNotUsed(request), vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get access to information and data
    auto in_info = input_vector[0]->GetInformationObject(0);
    auto input = vtkRectilinearGrid::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(input, "Input is not a rectilinear grid.");

    // Get output
    auto out_info = output_vector->GetInformationObject(0);
    auto output = vtkStructuredGrid::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(output, "Output is not a structured grid.");

    // Create structured grid from rectilinear grid
    std::array<int, 6> extent{};

    input->GetExtent(extent.data());

    auto coordinates_x = input->GetXCoordinates();
    auto coordinates_y = input->GetYCoordinates();
    auto coordinates_z = input->GetZCoordinates();

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(
        (static_cast<vtkIdType>(extent[1]) - extent[0] + 1) *
        (static_cast<vtkIdType>(extent[3]) - extent[2] + 1) *
        (static_cast<vtkIdType>(extent[5]) - extent[4] + 1));

    output->SetExtent(extent.data());

    vtkIdType index = 0;

    for (vtkIdType k = 0; k < coordinates_z->GetNumberOfTuples(); ++k)
    {
        for (vtkIdType j = 0; j < coordinates_y->GetNumberOfTuples(); ++j)
        {
            for (vtkIdType i = 0; i < coordinates_x->GetNumberOfTuples(); ++i)
            {
                const std::array<double, 3> point{
                    coordinates_x->GetComponent(i, 0),
                    coordinates_y->GetComponent(j, 0),
                    coordinates_z->GetComponent(k, 0) };

                points->SetPoint(index, point.data());

                ++index;
            }
        }
    }

    output->SetPoints(points);

    // Copy data
    output->GetPointData()->ShallowCopy(input->GetPointData());
    output->GetCellData()->ShallowCopy(input->GetCellData());
    output->GetFieldData()->ShallowCopy(input->GetFieldData());

    return 1;
}
