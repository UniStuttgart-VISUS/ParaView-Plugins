#include "regular_grid.h"

#include "vtkDoubleArray.h"
#include "vtkObjectFactory.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <stdexcept>

vtkStandardNewMacro(regular_grid);

regular_grid::regular_grid()
{
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(1);
}

int regular_grid::ProcessRequest(vtkInformation* request, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Create an output object of the correct type.
    if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA_OBJECT()))
    {
        return this->RequestDataObject(request, input_vector, output_vector);
    }

    // Generate the data
    if (request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
    {
        return this->RequestInformation(request, input_vector, output_vector);
    }

    if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
    {
        return this->RequestData(request, input_vector, output_vector);
    }

    if (request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
    {
        return this->RequestUpdateExtent(request, input_vector, output_vector);
    }

    return this->Superclass::ProcessRequest(request, input_vector, output_vector);
}

int regular_grid::RequestDataObject(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get VTK output type
    const auto output_type = (this->OutputType == 0) ? VTK_IMAGE_DATA : VTK_RECTILINEAR_GRID;

    // Create appropriate data object if the type differs from the previous execution
    auto output = output_vector->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT());

    if (output == nullptr || output->GetDataObjectType() != output_type)
    {
        // Delete previous output
        if (output != nullptr)
        {
            output->Delete();
        }

        // Create output object of appropriate type
        switch (output_type)
        {
        case VTK_RECTILINEAR_GRID:
            output = vtkRectilinearGrid::New();
            break;
        case VTK_IMAGE_DATA:
        default:
            output = vtkImageData::New();
            break;
        }

        // Set created object
        output_vector->GetInformationObject(0)->Set(vtkDataObject::DATA_OBJECT(), output);
        this->GetOutputPortInformation(0)->Set(vtkDataObject::DATA_EXTENT_TYPE(), output->GetExtentType());
    }

    return 1;
}

int regular_grid::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    vtkInformation* output_info = output_vector->GetInformationObject(0);

    // Set extent
    const int extent[6] = {
        0, this->Dimension[0],
        0, this->Dimension[1],
        0, this->Dimension[2]
    };

    output_info->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);

    // Cannot produce sub extents
    output_info->Set(CAN_PRODUCE_SUB_EXTENT(), 0);

    return 1;
}

int regular_grid::FillOutputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
        return 1;
    }

    return 0;
}

int regular_grid::RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*)
{
    return 1;
}

int regular_grid::RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    // Get output
    auto* out_info = output_vector->GetInformationObject(0);
    auto* output = out_info->Get(vtkDataObject::DATA_OBJECT());

    // Check parameters
    if (this->Dimension[0] <= 0 || this->Dimension[1] <= 0 || this->Dimension[2] <= 0)
    {
        std::cerr << "The number of cells must be larger than zero for each direction" << std::endl;
        return 0;
    }

    if (this->CellSize[0] <= 0.0 || this->CellSize[1] <= 0.0 || this->CellSize[2] <= 0.0)
    {
        std::cerr << "All cell sizes must be larger than zero" << std::endl;
        return 0;
    }

    // Calculate extent
    int extent[6] = {
        0, this->Dimension[0],
        0, this->Dimension[1],
        0, this->Dimension[2]
    };

    // Get origin depending on anchor
    const std::array<double, 3> origin{
        this->Origin[0] - ((this->Anchor == 0) ? 0.0 : ((this->CellSize[0] * this->Dimension[0]) / 2.0)),
        this->Origin[1] - ((this->Anchor == 0) ? 0.0 : ((this->CellSize[1] * this->Dimension[1]) / 2.0)),
        this->Origin[2] - ((this->Anchor == 0) ? 0.0 : ((this->CellSize[2] * this->Dimension[2]) / 2.0))
    };

    // Create regular grid
    if (this->OutputType == 0)
    {
        auto regular_grid = vtkImageData::SafeDownCast(output);

        if (regular_grid == nullptr)
        {
            std::cerr << "Data object type and requested type do not match." << std::endl;
            return 0;
        }

        // Set grid information
        regular_grid->SetExtent(extent);
        regular_grid->SetOrigin(origin.data());
        regular_grid->SetSpacing(this->CellSize.data());
    }

    // Create rectilinear grid
    if (this->OutputType == 1)
    {
        auto rectilinear_grid = vtkRectilinearGrid::SafeDownCast(output);

        if (rectilinear_grid == nullptr)
        {
            std::cerr << "Data object type and requested type do not match." << std::endl;
            return 0;
        }

        // Create arrays for node coordinates
        auto coordinates_x = vtkSmartPointer<vtkDoubleArray>::New();
        auto coordinates_y = vtkSmartPointer<vtkDoubleArray>::New();
        auto coordinates_z = vtkSmartPointer<vtkDoubleArray>::New();

        coordinates_x->SetNumberOfComponents(1);
        coordinates_y->SetNumberOfComponents(1);
        coordinates_z->SetNumberOfComponents(1);

        coordinates_x->SetNumberOfTuples(this->Dimension[0] + 1LL);
        coordinates_y->SetNumberOfTuples(this->Dimension[1] + 1LL);
        coordinates_z->SetNumberOfTuples(this->Dimension[2] + 1LL);

        for (vtkIdType i = 0; i <= this->Dimension[0]; ++i)
        {
            coordinates_x->SetValue(i, origin[0] + i * this->CellSize[0]);
        }
        for (vtkIdType i = 0; i <= this->Dimension[1]; ++i)
        {
            coordinates_y->SetValue(i, origin[1] + i * this->CellSize[1]);
        }
        for (vtkIdType i = 0; i <= this->Dimension[2]; ++i)
        {
            coordinates_z->SetValue(i, origin[2] + i * this->CellSize[2]);
        }

        // Set grid information
        rectilinear_grid->SetExtent(extent);
        rectilinear_grid->SetXCoordinates(coordinates_x);
        rectilinear_grid->SetYCoordinates(coordinates_y);
        rectilinear_grid->SetZCoordinates(coordinates_z);
    }

    return 1;
}
