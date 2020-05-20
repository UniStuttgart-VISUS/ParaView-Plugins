#include "filter_lines.h"

#include "common/math.h"

#include "vtkCellArray.h"
#include "vtkDataObject.h"
#include "vtkDoubleArray.h"
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

#include "Eigen/Dense"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

vtkStandardNewMacro(filter_lines);

filter_lines::filter_lines()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

int filter_lines::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
        return 1;
    }

    return 0;
}

int filter_lines::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    return 1;
}

int filter_lines::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get input lines
    auto* in_info = input_vector[0]->GetInformationObject(0);
    auto* input = vtkPolyData::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    std::vector<std::vector<vtkIdType>> lines(input->GetLines()->GetNumberOfCells());

    input->GetLines()->InitTraversal();
    for (auto& line : lines)
    {
        auto point_ids = vtkSmartPointer<vtkIdList>::New();
        input->GetLines()->GetNextCell(point_ids);

        line.resize(point_ids->GetNumberOfIds());

        for (std::size_t p = 0; p < line.size(); ++p)
        {
            line[p] = point_ids->GetId(p);
        }
    }

    // Check parameters
    if (this->MaxAngle < 0.0 || this->MaxAngle > 180.0)
    {
        std::cerr << "The maximum angle must not be negative or larger than 180 degrees" << std::endl;
        return 0;
    }

    if (this->MinSize < 0)
    {
        std::cerr << "The minimum size must not be negative" << std::endl;
        return 0;
    }

    // Filter by angle
    if (this->AngleFilter)
    {
        auto vectors = this->GetInputArrayToProcess(0, &input_vector[0]);

        if (vectors == nullptr)
        {
            std::cerr << "Cannot filter by angle: Not all arrays provided." << std::endl;
            return 0;
        }

        if (vectors->GetNumberOfTuples() != input->GetNumberOfPoints())
        {
            std::cerr << "Cannot filter by angle: Array must be provided as point data." << std::endl;
            return 0;
        }

        if (vectors->GetNumberOfComponents() != 3)
        {
            std::cerr << "Cannot filter by angle: Vectors must be three-dimensional." << std::endl;
            return 0;
        }

        const auto angle = pi * this->MaxAngle / 180.0;

        std::vector<std::vector<vtkIdType>> filtered_lines;
        filtered_lines.reserve(lines.size());

        for (const auto& line : lines)
        {
            filtered_lines.emplace_back();
            filtered_lines.back().reserve(line.size());

            // Look at each segment and calculate the angle
            for (std::size_t p = 0; p < line.size() - 1; ++p)
            {
                Eigen::Vector3d origin, target, vector;

                input->GetPoint(line[p], origin.data());
                input->GetPoint(line[p + 1], target.data());
                vectors->GetTuple(line[p], vector.data());

                const auto angle_1 = std::acos(vector.normalized().dot((target - origin).normalized()));
                const auto angle_2 = std::acos(vector.normalized().dot((origin - target).normalized()));

                if (std::min(angle_1, angle_2) < angle)
                {
                    // Add point to line
                    filtered_lines.back().push_back(line[p]);
                }
                else if (!filtered_lines.back().empty())
                {
                    // Start next line
                    filtered_lines.emplace_back();
                    filtered_lines.back().reserve(line.size() - p);
                }
            }

            // Always just add the last point
            filtered_lines.back().push_back(line.back());
        }

        std::swap(lines, filtered_lines);
    }

    // Filter lines by size and remove "lines" of one element
    std::size_t num_indices = 0;

    {
        std::vector<std::vector<vtkIdType>> filtered_lines;
        filtered_lines.reserve(lines.size());

        for (const auto& line : lines)
        {
            if (line.size() > (this->SizeFilter ? std::max(1, this->MinSize) : 1))
            {
                filtered_lines.push_back(line);
                num_indices += line.size();
            }
        }

        std::swap(lines, filtered_lines);
    }

    // Create output
    auto* out_info = output_vector->GetInformationObject(0);
    auto* output = vtkPolyData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    output->SetPoints(input->GetPoints());

    auto cell_indices = vtkSmartPointer<vtkIdTypeArray>::New();
    cell_indices->SetNumberOfComponents(1);
    cell_indices->Allocate(lines.size() + num_indices);

    for (const auto& line : lines)
    {
        cell_indices->InsertNextValue(line.size());

        for (const auto& index : line)
        {
            cell_indices->InsertNextValue(index);
        }
    }

    auto cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetCells(lines.size(), cell_indices);

    output->SetLines(cells);
    output->GetPointData()->ShallowCopy(input->GetPointData());

    return 1;
}
