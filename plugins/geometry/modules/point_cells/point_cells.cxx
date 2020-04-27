#include "point_cells.h"

#include "vtkCellArray.h"
#include "vtkDataObject.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

vtkStandardNewMacro(point_cells);

point_cells::point_cells()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

int point_cells::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
        return 1;
    }

    return 0;
}

int point_cells::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    return 1;
}

int point_cells::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get input and output
    auto* in_info = input_vector[0]->GetInformationObject(0);
    auto* input = vtkPointSet::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    auto* out_info = output_vector->GetInformationObject(0);
    auto* output = vtkPolyData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    // Copy points
    output->SetPoints(input->GetPoints());
    output->GetPointData()->ShallowCopy(input->GetPointData());

    // Create cells
    auto cell_indices = vtkSmartPointer<vtkIdTypeArray>::New();
    cell_indices->SetNumberOfComponents(1);
    cell_indices->Allocate(input->GetNumberOfPoints() * 2);

    for (vtkIdType index = 0; index < input->GetNumberOfPoints(); ++index)
    {
        cell_indices->InsertNextValue(1);
        cell_indices->InsertNextValue(index);
    }

    auto cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetCells(input->GetNumberOfPoints(), cell_indices);

    output->SetVerts(cells);

    return 1;
}
