#include "structured_grid_to_rectilinear_grid.h"

#include "common/checks.h"

#include "vtkCellData.h"
#include "vtkDoubleArray.h"
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

vtkStandardNewMacro(structured_grid_to_rectilinear_grid);

structured_grid_to_rectilinear_grid::structured_grid_to_rectilinear_grid()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

structured_grid_to_rectilinear_grid::~structured_grid_to_rectilinear_grid() {}

int structured_grid_to_rectilinear_grid::FillInputPortInformation(int port, vtkInformation* info)
{
    if (!this->Superclass::FillInputPortInformation(port, info))
    {
        return 0;
    }

    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkStructuredGrid");
        return 1;
    }

    return 0;
}

int structured_grid_to_rectilinear_grid::RequestUpdateExtent(vtkInformation* vtkNotUsed(request), vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    return 1;
}

int structured_grid_to_rectilinear_grid::RequestData(vtkInformation* vtkNotUsed(request), vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get access to information and data
    auto in_info = input_vector[0]->GetInformationObject(0);
    auto input = vtkStructuredGrid::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(input, "Input is not a structured grid.");

    // Get output
    auto out_info = output_vector->GetInformationObject(0);
    auto output = vtkRectilinearGrid::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(output, "Output is not a rectilinear grid.");

    // Create rectilinear grid from structured grid
    std::array<int, 6> extent{};

    input->GetExtent(extent.data());
    output->SetExtent(extent.data());

    auto points = input->GetPoints();

    auto coordinates_x = vtkSmartPointer<vtkDoubleArray>::New();
    auto coordinates_y = vtkSmartPointer<vtkDoubleArray>::New();
    auto coordinates_z = vtkSmartPointer<vtkDoubleArray>::New();

    coordinates_x->SetNumberOfComponents(1);
    coordinates_y->SetNumberOfComponents(1);
    coordinates_z->SetNumberOfComponents(1);

    coordinates_x->SetNumberOfTuples(static_cast<vtkIdType>(extent[1]) - extent[0] + 1);
    coordinates_y->SetNumberOfTuples(static_cast<vtkIdType>(extent[3]) - extent[2] + 1);
    coordinates_z->SetNumberOfTuples(static_cast<vtkIdType>(extent[5]) - extent[4] + 1);

    vtkIdType index = 0;

    for (vtkIdType k = extent[4]; k <= extent[5]; ++k)
    {
        for (vtkIdType j = extent[2]; j <= extent[3]; ++j)
        {
            for (vtkIdType i = extent[0]; i <= extent[1]; ++i)
            {
                std::array<double, 3> point{};

                points->GetPoint(index, point.data());

                if (j == 0 && k == 0)
                {
                    coordinates_x->SetValue(i, point[0]);
                }
                else if (std::abs(coordinates_x->GetValue(i) - point[0]) > 1.0e-6)
                {
                    std::cerr << "Input cannot be converted to rectilinear grid: X coordinates are not equal for every variation of Y, Z." << std::endl;

                    return 0;
                }

                if (i == 0 && k == 0)
                {
                    coordinates_y->SetValue(j, point[1]);
                }
                else if (std::abs(coordinates_y->GetValue(j) - point[1]) > 1.0e-6)
                {
                    std::cerr << "Input cannot be converted to rectilinear grid: Y coordinates are not equal for every variation of X, Z." << std::endl;

                    return 0;
                }

                if (i == 0 && j == 0)
                {
                    coordinates_z->SetValue(k, point[2]);
                }
                else if (std::abs(coordinates_z->GetValue(k) - point[2]) > 1.0e-6)
                {
                    std::cerr << "Input cannot be converted to rectilinear grid: Z coordinates are not equal for every variation of X, Y." << std::endl;

                    return 0;
                }

                ++index;
            }
        }
    }

    output->SetXCoordinates(coordinates_x);
    output->SetYCoordinates(coordinates_y);
    output->SetZCoordinates(coordinates_z);

    // Copy data
    output->GetPointData()->ShallowCopy(input->GetPointData());
    output->GetCellData()->ShallowCopy(input->GetCellData());
    output->GetFieldData()->ShallowCopy(input->GetFieldData());

    return 1;
}
