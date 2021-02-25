#include "grid_to_polydata.h"

#include "common/checks.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkImageData.h"
#include "vtkSmartPointer.h"
#include "vtkVertex.h"
#include "vtkLine.h"
#include "vtkPolygon.h"

#include <array>

vtkStandardNewMacro(grid_to_polydata);

grid_to_polydata::grid_to_polydata()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

grid_to_polydata::~grid_to_polydata() {}

int grid_to_polydata::FillInputPortInformation(int port, vtkInformation* info)
{
    if (!this->Superclass::FillInputPortInformation(port, info))
    {
        return 0;
    }

    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
        return 1;
    }

    return 0;
}

int grid_to_polydata::RequestUpdateExtent(vtkInformation* vtkNotUsed(request), vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    return 1;
}

int grid_to_polydata::RequestData(vtkInformation* vtkNotUsed(request), vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get access to information and data
    auto in_info = input_vector[0]->GetInformationObject(0);
    auto input = vtkUnstructuredGrid::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(input, "Input is not an unstructured grid.");

    // Get output
    auto out_info = output_vector->GetInformationObject(0);
    auto output = vtkPolyData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(output, "Output is not a polydata object.");

    // Copy points
    auto new_points = vtkSmartPointer<vtkPoints>::New();
    new_points->ShallowCopy(input->GetPoints());

    output->SetPoints(new_points);

    // Copy cells
    input->BuildLinks();

    auto cells = input->GetCells();
    auto cell_indices = cells->GetData();

    auto verts = vtkSmartPointer<vtkCellArray>::New();
    auto lines = vtkSmartPointer<vtkCellArray>::New();
    auto polys = vtkSmartPointer<vtkCellArray>::New();

    enum class state_t { expect_num, expect_val } state = state_t::expect_num;

    vtkIdType num_elements, num_element_index;
    vtkSmartPointer<vtkCell> vert, line, poly;
    int cell_type;

    for (vtkIdType index = 0, cell_index = 0; index < cell_indices->GetNumberOfValues(); ++index)
    {
        switch (state)
        {
        case state_t::expect_num:
            num_elements = cell_indices->GetValue(index);
            num_element_index = 0;

            cell_type = input->GetCellType(cell_index++);

            switch (cell_type)
            {
            case VTK_VERTEX:
                vert = vtkSmartPointer<vtkVertex>::New();
                vert->GetPointIds()->SetNumberOfIds(1);

                break;
            case VTK_LINE:
                line = vtkSmartPointer<vtkLine>::New();
                line->GetPointIds()->SetNumberOfIds(num_elements);

                break;
            case VTK_TRIANGLE:
            case VTK_POLYGON:
                poly = vtkSmartPointer<vtkPolygon>::New();
                poly->GetPointIds()->SetNumberOfIds(num_elements);
            }

            state = state_t::expect_val;

            break;
        case state_t::expect_val:
            switch (cell_type)
            {
            case VTK_VERTEX:
                vert->GetPointIds()->SetId(num_element_index, cell_indices->GetValue(index));

                break;
            case VTK_LINE:
                line->GetPointIds()->SetId(num_element_index, cell_indices->GetValue(index));

                break;
            case VTK_TRIANGLE:
            case VTK_POLYGON:
                poly->GetPointIds()->SetId(num_element_index, cell_indices->GetValue(index));
            }

            if (++num_element_index == num_elements)
            {
                switch (cell_type)
                {
                case VTK_VERTEX:
                    verts->InsertNextCell(vert);

                    break;
                case VTK_LINE:
                    lines->InsertNextCell(line);

                    break;
                case VTK_TRIANGLE:
                case VTK_POLYGON:
                    polys->InsertNextCell(poly);
                }

                state = state_t::expect_num;
            }

            break;
        }
    }

    output->SetVerts(verts);
    output->SetLines(lines);
    output->SetPolys(polys);

    // Copy arrays
    output->GetPointData()->ShallowCopy(input->GetPointData());
    output->GetCellData()->ShallowCopy(input->GetCellData());
    output->GetFieldData()->ShallowCopy(input->GetFieldData());

    return 1;
}
