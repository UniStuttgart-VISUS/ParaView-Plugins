#include "logarithmic_spiral.h"

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

vtkStandardNewMacro(logarithmic_spiral);

logarithmic_spiral::logarithmic_spiral()
{
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(1);
}

int logarithmic_spiral::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    return 1;
}

int logarithmic_spiral::RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    // Get output
    auto *out_info = output_vector->GetInformationObject(0);
    auto *output = vtkPolyData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    // Create parameterized "virtual" line starting at the origin, around which the logarithmic spiral is created
    const auto parameter_increment = this->Length / (this->NumPoints - 1);

    // Depending on the parameter, calculate position on the logarithmic spiral and set new point
    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(this->NumPoints);

    auto cell_ids = vtkSmartPointer<vtkIdTypeArray>::New();
    cell_ids->SetNumberOfValues(this->NumPoints + 1);
    cell_ids->SetValue(0, this->NumPoints);

    for (int i = 0; i < this->NumPoints; ++i)
    {
        std::array<double, 3> logarithmic_spiral_point;

        const auto spiral = this->SizeFactor * std::exp(this->Slope * i * parameter_increment);

        logarithmic_spiral_point[0] = spiral * std::cos(i * parameter_increment);
        logarithmic_spiral_point[1] = spiral * std::sin(i * parameter_increment);
        logarithmic_spiral_point[2] = this->Lift * i * parameter_increment;

        points->SetPoint(i, logarithmic_spiral_point.data());
        cell_ids->SetValue(i + 1, i);
    }

    // Create line
    auto cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetCells(1, cell_ids);

    output->SetPoints(points);
    output->SetLines(cells);

    return 1;
}
