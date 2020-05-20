#include "trigonometric_function.h"

#include "vtkCellArray.h"
#include "vtkDataObject.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <array>
#include <cmath>
#include <iostream>

vtkStandardNewMacro(trigonometric_function);

trigonometric_function::trigonometric_function()
{
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(1);
}

int trigonometric_function::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    return 1;
}

int trigonometric_function::RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    // Get output
    auto *out_info = output_vector->GetInformationObject(0);
    auto *output = vtkPolyData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    // Check parameters
    if (this->NumPoints <= 0)
    {
        std::cerr << "The number of output points must be larger than zero" << std::endl;
        return 0;
    }

    if (this->Length <= 0)
    {
        std::cerr << "The length must be larger than zero" << std::endl;
        return 0;
    }

    // Create parameterized "virtual" line starting at the origin, around which the trigonometric function is created
    const auto x_increment = this->Length / (this->NumPoints - 1);

    // Depending on the parameter in [0,1], calculate position on the trigonometric function and set new point
    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(this->NumPoints);

    auto cell_ids = vtkSmartPointer<vtkIdTypeArray>::New();
    cell_ids->SetNumberOfValues(this->NumPoints + 1);
    cell_ids->SetValue(0, this->NumPoints);

    for (int i = 0; i < this->NumPoints; ++i)
    {
        std::array<double, 3> trigonometric_function_point;
        trigonometric_function_point[0] = i * x_increment;
        trigonometric_function_point[2] = 0.0;

        if (this->Function == 0) // sine
        {
            trigonometric_function_point[1] = this->Alpha * std::sin(this->Beta * i * x_increment + this->Gamma) + this->Delta;
        }
        else if (this->Function == 1) // cosine
        {
            trigonometric_function_point[1] = this->Alpha * std::cos(this->Beta * i * x_increment + this->Gamma) + this->Delta;
        }
        else // tangens
        {
            trigonometric_function_point[1] = this->Alpha * std::tan(this->Beta * i * x_increment + this->Gamma) + this->Delta;
        }

        points->SetPoint(i, trigonometric_function_point.data());
        cell_ids->SetValue(i + 1, i);
    }

    // Create line
    auto cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetCells(1, cell_ids);

    output->SetPoints(points);
    output->SetLines(cells);

    return 1;
}
