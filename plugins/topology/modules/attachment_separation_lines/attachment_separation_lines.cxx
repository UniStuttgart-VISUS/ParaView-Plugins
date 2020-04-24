#include "attachment_separation_lines.h"

#include "vtkCellArray.h"
#include "vtkDataObject.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "Eigen/Dense"

#include <iostream>

vtkStandardNewMacro(attachment_separation_lines);

attachment_separation_lines::attachment_separation_lines()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

int attachment_separation_lines::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
        return 1;
    }

    return 0;
}

int attachment_separation_lines::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    return 1;
}

int attachment_separation_lines::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get input surface mesh
    auto* in_info = input_vector[0]->GetInformationObject(0);
    auto* input = vtkPolyData::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    // Get input arrays
    auto normals = this->GetInputArrayToProcess(0, &input_vector[0]);
    auto velocities = this->GetInputArrayToProcess(1, &input_vector[0]);
    auto Jacobians = this->GetInputArrayToProcess(2, &input_vector[0]);

    if (normals == nullptr ||
        velocities == nullptr ||
        Jacobians == nullptr)
    {
        std::cerr << "Not all arrays provided." << std::endl;
        return 0;
    }

    if (normals->GetNumberOfTuples() != input->GetNumberOfPoints() ||
        velocities->GetNumberOfTuples() != input->GetNumberOfPoints() ||
        Jacobians->GetNumberOfTuples() != input->GetNumberOfPoints())
    {
        std::cerr << "Arrays must be provided as point data." << std::endl;
        return 0;
    }

    if (normals->GetNumberOfComponents() != 3 ||
        velocities->GetNumberOfComponents() != 3 ||
        Jacobians->GetNumberOfComponents() != 9)
    {
        std::cerr << "Input must be three-dimensional." << std::endl;
        return 0;
    }

    // Compute
    auto values = vtkSmartPointer<vtkDoubleArray>::New();
    values->SetNumberOfComponents(1);
    values->SetNumberOfTuples(velocities->GetNumberOfTuples());
    values->SetName("Values");

    auto attachment_values = vtkSmartPointer<vtkDoubleArray>::New();
    attachment_values->SetNumberOfComponents(1);
    attachment_values->SetNumberOfTuples(velocities->GetNumberOfTuples());
    attachment_values->SetName("Attachment values");

    auto separation_values = vtkSmartPointer<vtkDoubleArray>::New();
    separation_values->SetNumberOfComponents(1);
    separation_values->SetNumberOfTuples(velocities->GetNumberOfTuples());
    separation_values->SetName("Separation values");

    auto attachment_eigenvectors = vtkSmartPointer<vtkDoubleArray>::New();
    attachment_eigenvectors->SetNumberOfComponents(3);
    attachment_eigenvectors->SetNumberOfTuples(velocities->GetNumberOfTuples());
    attachment_eigenvectors->SetName("Attachment eigenvectors");

    auto separation_eigenvectors = vtkSmartPointer<vtkDoubleArray>::New();
    separation_eigenvectors->SetNumberOfComponents(3);
    separation_eigenvectors->SetNumberOfTuples(velocities->GetNumberOfTuples());
    separation_eigenvectors->SetName("Separation eigenvectors");

    Eigen::Vector3d normal, velocity;
    Eigen::Matrix3d Jacobian;

    for (vtkIdType i = 0; i < velocities->GetNumberOfTuples(); ++i)
    {
        normals->GetTuple(i, normal.data());
        velocities->GetTuple(i, velocity.data());
        Jacobians->GetTuple(i, Jacobian.data());

        // Create orthonormal system on the surface
        Eigen::Vector3d basis_1, basis_2;

        if (normal[0] == 0.0)
        {
            basis_1[0] = 0.0;
            basis_1[1] = -normal[2];
            basis_1[2] = normal[1];
        }
        else if (normal[1] == 0.0)
        {
            basis_1[0] = -normal[2];
            basis_1[1] = 0.0;
            basis_1[2] = normal[0];
        }
        else
        {
            basis_1[0] = -normal[1];
            basis_1[1] = normal[0];
            basis_1[2] = 0.0;
        }

        basis_1.normalize();
        basis_2 = normal.cross(basis_1).normalized();

        // Project into the 2D system
        const Eigen::Vector2d projected_velocity(velocity.dot(basis_1), velocity.dot(basis_2));

        const auto w1 = Jacobian * basis_1;
        const auto w2 = Jacobian * basis_2;

        Eigen::Matrix2d projected_Jacobian;
        projected_Jacobian << w1.dot(basis_1), w2.dot(basis_1), w1.dot(basis_2), w2.dot(basis_2);

        // Compute value from acceleration
        const auto Ju = projected_Jacobian * projected_velocity;

        Eigen::Matrix2d uJu;
        uJu <<
            projected_velocity[0], Ju[0],
            projected_velocity[1], Ju[1];

        values->SetValue(i, uJu.determinant());

        // Compute values from Eigen decomposition, separately for attachment and separation
        const Eigen::EigenSolver<Eigen::Matrix2d> eigen(projected_Jacobian, true);
        const auto eigenvalues = eigen.eigenvalues();
        const auto eigenvectors = eigen.eigenvectors();

        if (eigenvalues.imag().isZero() && eigenvectors.imag().isZero())
        {
            Eigen::Vector2d larger, smaller;

            if (eigenvalues[0].real() > eigenvalues[1].real())
            {
                larger = eigenvectors.col(0).real();
                smaller = eigenvectors.col(1).real();
            }
            else
            {
                larger = eigenvectors.col(1).real();
                smaller = eigenvectors.col(0).real();
            }

            Eigen::Matrix2d ue_attachment;
            ue_attachment <<
                projected_velocity[0], smaller[0],
                projected_velocity[1], smaller[1];

            Eigen::Matrix2d ue_separation;
            ue_separation <<
                projected_velocity[0], larger[0],
                projected_velocity[1], larger[1];

            attachment_values->SetValue(i, ue_attachment.determinant());
            separation_values->SetValue(i, ue_separation.determinant());

            // Transform eigenvectors back to 3D space and store them
            Eigen::Matrix3d basis_transformation;
            basis_transformation << basis_1, basis_2, normal;

            const Eigen::Vector3d attachment_vector = basis_transformation * Eigen::Vector3d(smaller[0], smaller[1], 0.0);
            const Eigen::Vector3d separation_vector = basis_transformation * Eigen::Vector3d(larger[0], larger[1], 0.0);

            attachment_eigenvectors->SetTuple(i, attachment_vector.data());
            separation_eigenvectors->SetTuple(i, separation_vector.data());
        }
        else
        {
            // Set recognizable dummy values
            attachment_values->SetValue(i, std::nan(""));
            separation_values->SetValue(i, std::nan(""));

            attachment_eigenvectors->SetTuple3(i, 0.0, 0.0, 0.0);
            separation_eigenvectors->SetTuple3(i, 0.0, 0.0, 0.0);
        }
    }

    // Create output
    auto* out_info = output_vector->GetInformationObject(0);
    auto* output = vtkPolyData::SafeDownCast(out_info->Get(vtkDataObject::DATA_OBJECT()));

    output->ShallowCopy(input);
    output->GetPointData()->AddArray(values);
    output->GetPointData()->AddArray(attachment_values);
    output->GetPointData()->AddArray(separation_values);
    output->GetPointData()->AddArray(attachment_eigenvectors);
    output->GetPointData()->AddArray(separation_eigenvectors);

    return 1;
}
