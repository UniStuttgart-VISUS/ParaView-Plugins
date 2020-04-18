#include "helix.h"

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

vtkStandardNewMacro(helix);

helix::helix()
{
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(1);
}

int helix::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    return 1;
}

int helix::RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    // Get output
    auto *out_info = output_vector->GetInformationObject(0);
    auto *output = vtkPolyData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    // Create parameterized "virtual" line starting at the origin, around which the helix is created
    const auto height_increment = this->Length / (this->NumPoints - 1);

    // Set angle increment to get the number of windings specified
    const auto angle_increment = this->Windings * (std::atan(1.0) * 8.0) / (this->NumPoints - 1);

    // Depending on the parameter in [0,1], calculate position on the helix and set new point
    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(this->NumPoints);

    auto cell_ids = vtkSmartPointer<vtkIdTypeArray>::New();
    cell_ids->SetNumberOfValues(this->NumPoints + 1);
    cell_ids->SetValue(0, this->NumPoints);

    for (int i = 0; i < this->NumPoints; ++i)
    {
        const auto theta = i * angle_increment;

        const std::array<double, 3> helix_point{
            this->Radius * std::cos(theta),
            i * height_increment,
            this->Radius * std::sin(theta) };

        points->SetPoint(i, helix_point.data());
        cell_ids->SetValue(i + 1, i);
    }

    // Create line
    auto cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetCells(1, cell_ids);

    output->SetPoints(points);
    output->SetLines(cells);

    return 1;
}
