#include "sort_line_points.h"

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

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <vector>

vtkStandardNewMacro(sort_line_points);

sort_line_points::sort_line_points()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

int sort_line_points::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
        return 1;
    }

    return 0;
}

int sort_line_points::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    return 1;
}

int sort_line_points::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get input
    auto* in_info = input_vector[0]->GetInformationObject(0);
    auto* input = vtkPolyData::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    // Get output
    auto* out_info = output_vector->GetInformationObject(0);
    auto* output = vtkPolyData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    // Get lines and sort their points
    input->GetLines()->InitTraversal();
    auto point_ids = vtkSmartPointer<vtkIdList>::New();

    std::map<vtkIdType, vtkIdType> point_map;
    vtkIdType new_index = 0;

    while (input->GetLines()->GetNextCell(point_ids))
    {
        // Create point map
        for (vtkIdType i = 0; i < point_ids->GetNumberOfIds(); ++i)
        {
            if (point_map.find(point_ids->GetId(i)) == point_map.end())
            {
                point_map[point_ids->GetId(i)] = new_index++;
            }
        }
    }

    // Copy points to new position
    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(point_map.size());

    for (const auto indices : point_map)
    {
        const auto old_index = indices.first;
        const auto new_index = indices.second;

        std::array<double, 3> point;
        input->GetPoint(old_index, point.data());
        points->SetPoint(new_index, point.data());
    }

    output->SetPoints(points);

    // Create cells
    auto cell_indices = input->GetLines()->GetData();
    auto cell_index = 0;
    auto num_cells = 0;

    while (cell_index < cell_indices->GetNumberOfValues())
    {
        const vtkIdType num_points = cell_indices->GetValue(cell_index++);

        for (vtkIdType i = 0; i < num_points; ++i, ++cell_index)
        {
            cell_indices->SetValue(cell_index, point_map[cell_indices->GetValue(cell_index)]);
        }

        ++num_cells;
    }

    auto cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetCells(num_cells, cell_indices);

    output->SetLines(cells);

    return 1;
}
