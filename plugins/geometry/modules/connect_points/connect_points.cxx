#include "connect_points.h"

#include "vtkCellArray.h"
#include "vtkDataObject.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

vtkStandardNewMacro(connect_points);

connect_points::connect_points()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

int connect_points::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
        return 1;
    }

    return 0;
}

int connect_points::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    return 1;
}

int connect_points::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get input and output
    auto* in_info = input_vector[0]->GetInformationObject(0);
    auto* input = vtkPolyData::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    auto* out_info = output_vector->GetInformationObject(0);
    auto* output = vtkPolyData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    // Copy points
    output->SetPoints(input->GetPoints());
    output->GetPointData()->ShallowCopy(input->GetPointData());

    // Create cells
    auto cell_indices = vtkSmartPointer<vtkIdTypeArray>::New();
    cell_indices->SetNumberOfComponents(1);
    cell_indices->Allocate(input->GetNumberOfPoints() + 1);

    cell_indices->InsertNextValue(input->GetNumberOfPoints());

    for (vtkIdType index = 0; index < input->GetNumberOfPoints(); ++index)
    {
        cell_indices->InsertNextValue(index);
    }

    auto cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetCells(1, cell_indices);

    output->SetLines(cells);

    return 1;
}
