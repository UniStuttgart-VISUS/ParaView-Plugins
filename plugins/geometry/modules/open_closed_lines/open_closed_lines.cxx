#include "open_closed_lines.h"

#include "vtkCellData.h"
#include "vtkDataObject.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <array>

vtkStandardNewMacro(open_closed_lines);

open_closed_lines::open_closed_lines()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

int open_closed_lines::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
        return 1;
    }

    return 0;
}

int open_closed_lines::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    return 1;
}

int open_closed_lines::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get input
    auto* in_info = input_vector[0]->GetInformationObject(0);
    auto* input = vtkPolyData::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    // Get output
    auto* out_info = output_vector->GetInformationObject(0);
    auto* output = vtkPolyData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    // Open closed lines depending on method
    auto cell_ids = vtkSmartPointer<vtkIdTypeArray>::New();
    cell_ids->SetNumberOfComponents(1);
    cell_ids->Allocate(input->GetLines()->GetData()->GetNumberOfValues());

    if (this->Method == 0)
    {
        // Remove last segment
        input->GetLines()->InitTraversal();

        auto point_ids = vtkSmartPointer<vtkIdList>::New();

        while (input->GetLines()->GetNextCell(point_ids))
        {
            vtkIdType new_size = point_ids->GetNumberOfIds();

            if (point_ids->GetId(0) == point_ids->GetId(point_ids->GetNumberOfIds() - 1))
            {
                --new_size;
            }

            cell_ids->InsertNextValue(new_size);

            for (vtkIdType i = 0; i < new_size; ++i)
            {
                cell_ids->InsertNextValue(point_ids->GetId(i));
            }
        }

        // Copy points
        output->SetPoints(input->GetPoints());
    }
    else
    {
        // Duplicate shared point
        input->GetLines()->InitTraversal();

        auto point_ids = vtkSmartPointer<vtkIdList>::New();

        auto new_points = vtkSmartPointer<vtkPoints>::New();
        new_points->DeepCopy(input->GetPoints());

        while (input->GetLines()->GetNextCell(point_ids))
        {
            cell_ids->InsertNextValue(point_ids->GetNumberOfIds());

            for (vtkIdType i = 0; i < point_ids->GetNumberOfIds() - 1; ++i)
            {
                cell_ids->InsertNextValue(point_ids->GetId(i));
            }

            if (point_ids->GetId(0) == point_ids->GetId(point_ids->GetNumberOfIds() - 1))
            {
                std::array<double, 3> point;
                new_points->GetPoint(point_ids->GetId(0), point.data());

                cell_ids->InsertNextValue(new_points->InsertNextPoint(point.data()));
            }
            else
            {
                cell_ids->InsertNextValue(point_ids->GetId(point_ids->GetNumberOfIds() - 1));
            }
        }

        // Set points
        output->SetPoints(new_points);
    }

    // Create cells
    auto cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetCells(input->GetLines()->GetNumberOfCells(), cell_ids);

    output->SetLines(cells);

    return 1;
}
