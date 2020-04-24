#include "connect_lines.h"

#include "vtkCellArray.h"
#include "vtkDataObject.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <iostream>
#include <list>

vtkStandardNewMacro(connect_lines);

connect_lines::connect_lines()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

int connect_lines::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
        return 1;
    }

    return 0;
}

int connect_lines::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    return 1;
}

int connect_lines::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get input lines
    auto* in_info = input_vector[0]->GetInformationObject(0);
    auto* input = vtkPolyData::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    std::list<std::list<vtkIdType>> lines(input->GetLines()->GetNumberOfCells());

    input->GetLines()->InitTraversal();
    for (auto& line : lines)
    {
        auto point_ids = vtkSmartPointer<vtkIdList>::New();
        input->GetLines()->GetNextCell(point_ids);

        for (std::size_t p = 0; p < static_cast<std::size_t>(point_ids->GetNumberOfIds()); ++p)
        {
            line.push_back(point_ids->GetId(p));
        }
    }

    // Create polylines
    std::list<std::list<vtkIdType>> polylines;

    std::size_t num_indices = 0;

    while (!lines.empty())
    {
        polylines.push_back(lines.front());
        lines.pop_front();

        bool change = true;

        while (change)
        {
            change = false;

            for (auto it = lines.begin(); it != lines.end();)
            {
                if (polylines.back().front() == it->front())
                {
                    polylines.back().push_front(it->back());
                    lines.erase(it++);

                    change = true;
                }
                else if (polylines.back().front() == it->back())
                {
                    polylines.back().push_front(it->front());
                    lines.erase(it++);

                    change = true;
                }
                else if (polylines.back().back() == it->front())
                {
                    polylines.back().push_back(it->back());
                    lines.erase(it++);

                    change = true;
                }
                else if (polylines.back().back() == it->back())
                {
                    polylines.back().push_back(it->front());
                    lines.erase(it++);

                    change = true;
                }
                else
                {
                    ++it;
                }
            }
        }

        num_indices += polylines.back().size();
    }

    // Create output
    auto* out_info = output_vector->GetInformationObject(0);
    auto* output = vtkPolyData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    output->SetPoints(input->GetPoints());

    auto cell_indices = vtkSmartPointer<vtkIdTypeArray>::New();
    cell_indices->SetNumberOfComponents(1);
    cell_indices->Allocate(polylines.size() + num_indices);

    for (const auto& line : polylines)
    {
        cell_indices->InsertNextValue(line.size());

        for (const auto& index : line)
        {
            cell_indices->InsertNextValue(index);
        }
    }

    auto cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetCells(polylines.size(), cell_indices);

    output->SetLines(cells);
    output->GetPointData()->ShallowCopy(input->GetPointData());

    return 1;
}
