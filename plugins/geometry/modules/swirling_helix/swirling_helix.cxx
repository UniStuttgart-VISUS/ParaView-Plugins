#include "swirling_helix.h"

#include "vtkCellArray.h"
#include "vtkDataObject.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "Eigen/Dense"

#include <cmath>

vtkStandardNewMacro(swirling_helix);

swirling_helix::swirling_helix()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

int swirling_helix::FillInputPortInformation(int port, vtkInformation* info)
{
    if (!this->Superclass::FillInputPortInformation(port, info))
    {
        return 0;
    }

    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
        return 1;
    }

    return 0;
}

int swirling_helix::RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*)
{
    return 1;
}

int swirling_helix::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get input
    auto* in_info = input_vector[0]->GetInformationObject(0);
    auto* input = vtkPolyData::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    auto* tangents = this->GetInputArrayToProcess(0, input_vector);
    auto* normals = this->GetInputArrayToProcess(1, input_vector);

    // Get output
    auto* out_info = output_vector->GetInformationObject(0);
    auto* output = vtkPolyData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    // Create output point and cell arrays
    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(input->GetPoints()->GetNumberOfPoints());

    auto cell_ids = vtkSmartPointer<vtkIdTypeArray>::New();
    cell_ids->Allocate(input->GetPoints()->GetNumberOfPoints() + input->GetLines()->GetNumberOfCells());

    // For each input line
    auto lines = input->GetLines();

    vtkIdType point_index = 0;

    for (vtkIdType l = 0; l < lines->GetNumberOfCells(); ++l)
    {
        auto point_list = vtkSmartPointer<vtkIdList>::New();
        lines->GetCell(l, point_list);

        cell_ids->InsertNextValue(point_list->GetNumberOfIds());

        // Set increment
        const auto angle_increment = this->Windings * (std::atan(1.0) * 8.0) / (point_list->GetNumberOfIds() - 1);

        // For each point
        for (vtkIdType p = 0; p < point_list->GetNumberOfIds(); ++p, ++point_index)
        {
            Eigen::Vector3d point;
            input->GetPoints()->GetPoint(point_list->GetId(p), point.data());

            // Create swirling helix point...
            Eigen::Vector3d tangent, normal;
            tangents->GetTuple(point_list->GetId(p), tangent.data());
            normals->GetTuple(point_list->GetId(p), normal.data());

            // ... create point relative to origin
            const Eigen::Vector3d point_relative_to_origin = this->Radius * normal;

            // ... rotate around the tangent
            const auto theta = p * angle_increment;

            const Eigen::Vector3d rotated_point_relative_to_origin =
                tangent.dot(point_relative_to_origin) * tangent
                + (std::cos(theta) * tangent.cross(point_relative_to_origin)).cross(tangent)
                + std::sin(theta) * tangent.cross(point_relative_to_origin);

            // ... add to original point
            point = point + rotated_point_relative_to_origin;

            // Set output point
            points->SetPoint(point_index, point.data());
            cell_ids->InsertNextValue(point_index);
        }
    }

    // Create cells
    auto cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetCells(lines->GetNumberOfCells(), cell_ids);

    output->SetPoints(points);
    output->SetLines(cells);

    return 1;
}
