#include "remove_long_segments.h"

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
#include <vector>

vtkStandardNewMacro(remove_long_segments);

remove_long_segments::remove_long_segments()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

int remove_long_segments::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
        return 1;
    }

    return 0;
}

int remove_long_segments::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    return 1;
}

int remove_long_segments::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get input
    auto* in_info = input_vector[0]->GetInformationObject(0);
    auto* input = vtkPolyData::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    // Get output
    auto* out_info = output_vector->GetInformationObject(0);
    auto* output = vtkPolyData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    // Get lines and remove long segments
    std::vector<std::vector<vtkIdType>> lines;

    input->GetLines()->InitTraversal();
    auto point_ids = vtkSmartPointer<vtkIdList>::New();

    while (input->GetLines()->GetNextCell(point_ids))
    {
        // Get median length
        std::vector<double> lengths(point_ids->GetNumberOfIds() - 1);

        for (vtkIdType i = 0; i < point_ids->GetNumberOfIds() - 1; ++i)
        {
            std::array<double, 3> point_1, point_2;
            input->GetPoint(point_ids->GetId(i), point_1.data());
            input->GetPoint(point_ids->GetId(i + 1), point_2.data());

            lengths[i] = std::sqrt(std::pow(point_1[0] - point_2[0], 2.0) + std::pow(point_1[1] - point_2[1], 2.0) + std::pow(point_1[2] - point_2[2], 2.0));
        }

        std::nth_element(lengths.begin(), lengths.begin() + lengths.size() / 2, lengths.end());
        const auto median_length = lengths[lengths.size() / 2];

        // Remove segments above the threshold
        const auto threshold = this->Absolute ? this->LengthFactor : (this->LengthFactor * median_length);

        // Remove segments above the threshold
        lines.emplace_back();
        lines.back().push_back(point_ids->GetId(0));

        for (vtkIdType i = 0; i < point_ids->GetNumberOfIds() - 1; ++i)
        {
            std::array<double, 3> point_1, point_2;
            input->GetPoint(point_ids->GetId(i), point_1.data());
            input->GetPoint(point_ids->GetId(i + 1), point_2.data());

            const auto length = std::sqrt(std::pow(point_1[0] - point_2[0], 2.0) + std::pow(point_1[1] - point_2[1], 2.0) + std::pow(point_1[2] - point_2[2], 2.0));

            if (length > threshold && !lines.back().empty())
            {
                lines.emplace_back();
            }

            lines.back().push_back(point_ids->GetId(i + 1));
        }
    }

    // Count number of cells and points
    std::size_t num_points = 0;
    std::size_t num_cells = 0;

    for (const auto& line : lines)
    {
        if (line.size() > 1)
        {
            ++num_cells;
            num_points += line.size();
        }
    }

    // Create cells
    auto cell_ids = vtkSmartPointer<vtkIdTypeArray>::New();
    cell_ids->SetNumberOfComponents(1);
    cell_ids->Allocate(static_cast<vtkIdType>(num_cells + num_points));

    for (const auto& line : lines)
    {
        if (line.size() > 1)
        {
            cell_ids->InsertNextValue(static_cast<vtkIdType>(line.size()));

            for (const auto point_id : line)
            {
                cell_ids->InsertNextValue(point_id);
            }
        }
    }

    auto cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetCells(static_cast<vtkIdType>(num_cells), cell_ids);

    output->SetLines(cells);

    // Copy points
    if (input->GetPoints() != nullptr)
    {
        auto points = vtkSmartPointer<vtkPoints>::New();
        points->ShallowCopy(input->GetPoints());

        output->SetPoints(points);
    }

    return 1;
}
