#include "smooth_lines.h"

#include "vtkDataObject.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "Eigen/Dense"

#include <vector>

vtkStandardNewMacro(smooth_lines);

smooth_lines::smooth_lines()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

int smooth_lines::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
        return 1;
    }

    return 0;
}

int smooth_lines::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    return 1;
}

int smooth_lines::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get input
    auto* in_info = input_vector[0]->GetInformationObject(0);
    auto* input = vtkPolyData::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    // Get output
    auto* out_info = output_vector->GetInformationObject(0);
    auto* output = vtkPolyData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    // Copy input to output and work on that copy
    output->DeepCopy(input);

    // Line-wise
    output->GetLines()->InitTraversal();

    auto point_indices = vtkSmartPointer<vtkIdList>::New();

    while (output->GetLines()->GetNextCell(point_indices))
    {
        std::vector<Eigen::Vector3d> points(point_indices->GetNumberOfIds()), temp_points(point_indices->GetNumberOfIds());

        for (vtkIdType index = 0; index < point_indices->GetNumberOfIds(); ++index)
        {
            output->GetPoints()->GetPoint(point_indices->GetId(index), points[index].data());
        }

        // Apply Taubin smoothing
        for (std::size_t i = 0; i < this->NumIterations; ++i)
        {
            for (std::size_t index = 0; index < points.size(); ++index)
            {
                temp_points[index] = gaussian_smoothing(points, index, this->Lambda);
            }

            for (std::size_t index = 0; index < points.size(); ++index)
            {
                points[index] = gaussian_smoothing(temp_points, index, this->Mu);
            }
        }

        // Write smoothed points to output
        for (vtkIdType index = 0; index < point_indices->GetNumberOfIds(); ++index)
        {
            output->GetPoints()->SetPoint(point_indices->GetId(index), points[index].data());
        }
    }

    return 1;
}

Eigen::Vector3d smooth_lines::gaussian_smoothing(const std::vector<Eigen::Vector3d>& points, const std::size_t index, const double weight) const
{
    auto point = points[index];

    if (index == 0)
    {
        point = point + weight * (points[index + 1] - point);
    }
    else if (index == points.size() - 1)
    {
        point = point + weight * (points[index - 1] - point);
    }
    else
    {
        point = point + weight * (
            0.5 * (points[index - 1] - point) +
            0.5 * (points[index + 1] - point));
    }

    return point;
}
