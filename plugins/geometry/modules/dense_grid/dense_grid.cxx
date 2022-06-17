#include "dense_grid.h"

#include "vtkCellArray.h"
#include "vtkDataObject.h"
#include "vtkDataSet.h"
#include "vtkIdTypeArray.h"
#include "vtkImageData.h"
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
#include <iostream>

vtkStandardNewMacro(dense_grid);

dense_grid::dense_grid()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

int dense_grid::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
        return 1;
    }

    return 0;
}

int dense_grid::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    return 1;
}

int dense_grid::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get input and output
    auto* in_info = input_vector[0]->GetInformationObject(0);
    auto* input = vtkDataSet::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    auto* out_info = output_vector->GetInformationObject(0);
    auto* output = vtkPolyData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    // Depending on input, get information
    auto* input_image = vtkImageData::SafeDownCast(input);
    auto* input_rectilinear_grid = vtkRectilinearGrid::SafeDownCast(input);
    auto* input_structured_grid = vtkStructuredGrid::SafeDownCast(input);

    std::array<int, 3> dimension{};

    if (input_image != nullptr)
    {
        input_image->GetDimensions(dimension.data());
    }
    else if (input_rectilinear_grid != nullptr)
    {
        input_rectilinear_grid->GetDimensions(dimension.data());
    }
    else if (input_structured_grid != nullptr)
    {
        input_structured_grid->GetDimensions(dimension.data());
    }

    const auto twoD = dimension[2] == 1;
    const auto oneD = dimension[1] == 1;

    // Create points
    std::array<double, 3> point{};

    if (input_image != nullptr)
    {
        std::array<double, 3> origin{}, spacing{};

        input_image->GetOrigin(origin.data());
        input_image->GetSpacing(spacing.data());

        auto points = vtkSmartPointer<vtkPoints>::New();

        for (std::size_t k = 0; k < dimension[2]; ++k)
        {
            for (std::size_t j = 0; j < dimension[1]; ++j)
            {
                for (std::size_t i = 0; i < dimension[0]; ++i)
                {
                    point[0] = origin[0] + i * spacing[0];
                    point[1] = origin[1] + j * spacing[1];
                    point[2] = origin[2] + k * spacing[2];

                    points->InsertNextPoint(point.data());
                }
            }
        }

        output->SetPoints(points);
    }
    else if (input_rectilinear_grid != nullptr)
    {
        auto x_coordinates = input_rectilinear_grid->GetXCoordinates();
        auto y_coordinates = input_rectilinear_grid->GetYCoordinates();
        auto z_coordinates = input_rectilinear_grid->GetZCoordinates();

        auto points = vtkSmartPointer<vtkPoints>::New();

        for (std::size_t k = 0; k < dimension[2]; ++k)
        {
            point[2] = z_coordinates->GetComponent(k, 0);

            for (std::size_t j = 0; j < dimension[1]; ++j)
            {
                point[1] = y_coordinates->GetComponent(j, 0);

                for (std::size_t i = 0; i < dimension[0]; ++i)
                {
                    point[0] = x_coordinates->GetComponent(i, 0);

                    points->InsertNextPoint(point.data());
                }
            }
        }

        output->SetPoints(points);
    }
    else if (input_structured_grid != nullptr)
    {
        auto points = vtkSmartPointer<vtkPoints>::New();

        points->ShallowCopy(input_structured_grid->GetPoints());
    }
    else
    {
        std::cerr << "Input data structure not supported." << std::endl;
        return 0;
    }

    // Create line segments
    auto lines = vtkSmartPointer<vtkCellArray>::New();
    std::array<vtkIdType, 2> line{};

    auto ratio = this->Ratio;
    if (ratio < 1)
    {
        std::cerr << "Warning: Sampling ratio must be at least 1." << std::endl;
        ratio = 1;
    }

    for (std::size_t k = 0; k < dimension[2];)
    {
        for (std::size_t j = 0; j < dimension[1];)
        {
            for (std::size_t i = 0; i < dimension[0];)
            {
                const vtkIdType index = i + dimension[0] * (j + dimension[1] * k);

                if (i < static_cast<long long>(dimension[0]) - 1)
                {
                    const vtkIdType index_right = std::min<long long>(i + ratio, static_cast<long long>(dimension[0]) - 1) + dimension[0] * (j + dimension[1] * k);

                    line = { index, index_right };
                    lines->InsertNextCell(2, line.data());
                }

                if (!oneD && j < static_cast<long long>(dimension[1]) - 1)
                {
                    const vtkIdType index_top = i + dimension[0] * (std::min<long long>(j + ratio, static_cast<long long>(dimension[1]) - 1) + dimension[1] * k);

                    line = { index, index_top };
                    lines->InsertNextCell(2, line.data());
                }

                if (!twoD && k < static_cast<long long>(dimension[2]) - 1)
                {
                    const vtkIdType index_back = i + dimension[0] * (j + dimension[1] * std::min<long long>(k + ratio, static_cast<long long>(dimension[2]) - 1));

                    line = { index, index_back };
                    lines->InsertNextCell(2, line.data());
                }

                if (i < static_cast<long long>(dimension[0]) - 1 && i + ratio >= dimension[0])
                {
                    i = static_cast<long long>(dimension[0]) - 1;
                }
                else
                {
                    i += ratio;
                }
            }

            if (j < static_cast<long long>(dimension[1]) - 1 && j + ratio >= dimension[1])
            {
                j = static_cast<long long>(dimension[1]) - 1;
            }
            else
            {
                j += ratio;
            }
        }

        if (k < static_cast<long long>(dimension[2]) - 1 && k + ratio >= dimension[2])
        {
            k = static_cast<long long>(dimension[2]) - 1;
        }
        else
        {
            k += ratio;
        }
    }

    output->SetLines(lines);

    output->GetPointData()->ShallowCopy(input->GetPointData());

    return 1;
}
