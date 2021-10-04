#include "structured_grid_to_unstructured_grid.h"

#include "common/checks.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"

#include <algorithm>
#include <array>

vtkStandardNewMacro(structured_grid_to_unstructured_grid);

structured_grid_to_unstructured_grid::structured_grid_to_unstructured_grid()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

structured_grid_to_unstructured_grid::~structured_grid_to_unstructured_grid() {}

int structured_grid_to_unstructured_grid::FillInputPortInformation(int port, vtkInformation* info)
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

int structured_grid_to_unstructured_grid::RequestUpdateExtent(vtkInformation* vtkNotUsed(request), vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    return 1;
}

int structured_grid_to_unstructured_grid::RequestData(vtkInformation* vtkNotUsed(request), vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get access to information and data
    auto in_info = input_vector[0]->GetInformationObject(0);
    auto input = vtkStructuredGrid::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(input, "Input is not a structured grid.");

    // Get output
    auto out_info = output_vector->GetInformationObject(0);
    auto output = vtkUnstructuredGrid::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(output, "Output is not an unstructured grid.");

    // Create unstructured grid from structured grid
    std::array<int, 6> extent{};

    input->GetExtent(extent.data());

    const auto dim_x = static_cast<vtkIdType>(extent[1]) - extent[0] + 1;
    const auto dim_y = static_cast<vtkIdType>(extent[3]) - extent[2] + 1;
    const auto dim_z = static_cast<vtkIdType>(extent[5]) - extent[4] + 1;

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->ShallowCopy(input->GetPoints());

    output->SetPoints(points);
    output->Allocate((dim_x - 1) * (dim_y - 1) * std::max(dim_z - 1, 1LL));

    auto get_index = [dim_x, dim_y, dim_z](const vtkIdType i, const vtkIdType j, const vtkIdType k) -> vtkIdType
    {
        return i + dim_x * (j + dim_y * k);
    };

    for (vtkIdType k = 0; k < std::max(dim_z - 1, 1LL); ++k)
    {
        for (vtkIdType j = 0; j < dim_y - 1; ++j)
        {
            for (vtkIdType i = 0; i < dim_x - 1; ++i)
            {
                const auto index = get_index(i, j, k);
                const auto index_right = get_index(i + 1, j, k);
                const auto index_top = get_index(i, j + 1, k);
                const auto index_right_top = get_index(i + 1, j + 1, k);

                if (dim_z == 1)
                {
                    const std::array<vtkIdType, 4> point_ids{ index, index_right, index_right_top, index_top };

                    output->InsertNextCell(VTK_QUAD, 4, point_ids.data());
                }
                else
                {
                    const auto index_back = get_index(i, j, k + 1);
                    const auto index_right_back = get_index(i + 1, j, k + 1);
                    const auto index_top_back = get_index(i, j + 1, k + 1);
                    const auto index_right_top_back = get_index(i + 1, j + 1, k + 1);

                    const std::array<vtkIdType, 8> point_ids{ index, index_right, index_right_top, index_top };

                    std::array<vtkIdType, 4> face0{ index, index_right, index_right_top, index_top }; // front
                    std::array<vtkIdType, 4> face1{ index_top_back, index_right_top_back, index_right_back, index_back }; // back
                    std::array<vtkIdType, 4> face2{ index_back, index_right_back, index_right, index }; // bottom
                    std::array<vtkIdType, 4> face3{ index_top, index_right_top, index_right_top_back, index_top_back }; // top
                    std::array<vtkIdType, 4> face4{ index, index_top, index_top_back, index_back }; // left
                    std::array<vtkIdType, 4> face5{ index_right, index_right_back, index_right_top_back, index_right_top }; // right

                    auto faces = vtkSmartPointer<vtkCellArray>::New();
                    faces->InsertNextCell(4, face0.data());
                    faces->InsertNextCell(4, face1.data());
                    faces->InsertNextCell(4, face2.data());
                    faces->InsertNextCell(4, face3.data());
                    faces->InsertNextCell(4, face4.data());
                    faces->InsertNextCell(4, face5.data());

                    output->InsertNextCell(VTK_POLYHEDRON, 8, point_ids.data(), 6, faces->GetData()->GetPointer(0));
                }
            }
        }
    }

    // Copy data
    output->GetPointData()->ShallowCopy(input->GetPointData());
    output->GetCellData()->ShallowCopy(input->GetCellData());
    output->GetFieldData()->ShallowCopy(input->GetFieldData());

    return 1;
}
