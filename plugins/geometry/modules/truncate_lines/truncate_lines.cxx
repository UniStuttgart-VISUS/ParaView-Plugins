#include "truncate_lines.h"

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
#include <unordered_set>
#include <vector>

vtkStandardNewMacro(truncate_lines);

truncate_lines::truncate_lines()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

int truncate_lines::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
        return 1;
    }

    return 0;
}

int truncate_lines::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    return 1;
}

int truncate_lines::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get input
    auto* in_info = input_vector[0]->GetInformationObject(0);
    auto* input = vtkPolyData::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    // Get output
    auto* out_info = output_vector->GetInformationObject(0);
    auto* output = vtkPolyData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    // Get lines and truncate them
    std::unordered_set<vtkIdType> point_indices;
    std::vector<std::vector<vtkIdType>> lines;

    std::size_t num_points = 0;

    if (this->NumPoints > 0)
    {
        input->GetLines()->InitTraversal();

        auto point_ids = vtkSmartPointer<vtkIdList>::New();

        while (input->GetLines()->GetNextCell(point_ids))
        {
            const auto size = point_ids->GetNumberOfIds();

            if (this->Offset < size)
            {
                lines.emplace_back();

                for (vtkIdType i = std::max(0, this->Offset); i < this->NumPoints && i < point_ids->GetNumberOfIds(); ++i)
                {
                    point_indices.insert(point_ids->GetId(i));
                    lines.back().push_back(point_ids->GetId(i));
                }

                num_points += lines.back().size();
            }
        }
    }

    // Create cells
    auto cell_ids = vtkSmartPointer<vtkIdTypeArray>::New();
    cell_ids->SetNumberOfComponents(1);
    cell_ids->Allocate(static_cast<vtkIdType>(lines.size() + num_points));

    for (const auto& line : lines)
    {
        cell_ids->InsertNextValue(static_cast<vtkIdType>(line.size()));

        for (const auto point_id : line)
        {
            cell_ids->InsertNextValue(point_id);
        }
    }

    auto cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetCells(static_cast<vtkIdType>(lines.size()), cell_ids);

    output->SetLines(cells);

    // Copy points
    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(point_indices.size());

    vtkIdType index = 0;

    for (const auto point_id : point_indices)
    {
        std::array<double, 3> point;
        input->GetPoint(point_id, point.data());

        points->SetPoint(index++, point.data());
    }

    output->SetPoints(points);

    return 1;
}
