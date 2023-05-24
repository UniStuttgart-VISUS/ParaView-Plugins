#include "separate_block.h"

#include "common/checks.h"

#include "vtkDataObject.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkRectilinearGrid.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStructuredGrid.h"
#include "vtkStructuredPoints.h"
#include "vtkType.h"
#include "vtkUnstructuredGrid.h"

#include <iostream>

vtkStandardNewMacro(separate_block);

separate_block::separate_block()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

int separate_block::ProcessRequest(vtkInformation* request, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
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

int separate_block::RequestDataObject(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get input block
    auto in_info = input_vector[0]->GetInformationObject(0);
    auto input = vtkMultiBlockDataSet::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(input, "Input is not a multiblock dataset.");

    if (this->BlockID < 0 || this->BlockID >= static_cast<int>(input->GetNumberOfBlocks()))
    {
        return 0;
    }

    // Create appropriate data object if the type differs from the previous execution
    auto output = output_vector->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT());

    if (output == nullptr || output->GetDataObjectType() != input->GetBlock(this->BlockID)->GetDataObjectType())
    {
        // Create output object of appropriate type
        auto* newOutput = input->GetBlock(this->BlockID)->NewInstance();
        output_vector->GetInformationObject(0)->Set(vtkDataObject::DATA_OBJECT(), newOutput);
        newOutput->Delete();
    }

    return 1;
}

int separate_block::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    return 1;
}

int separate_block::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
        return 1;
    }

    return 0;
}

int separate_block::FillOutputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
        return 1;
    }

    return 0;
}

int separate_block::RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*)
{
    return 1;
}

int separate_block::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get access to information and data
    auto in_info = input_vector[0]->GetInformationObject(0);
    auto input = vtkMultiBlockDataSet::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(input, "Input is not a multiblock dataset.");

    // Get output
    auto out_info = output_vector->GetInformationObject(0);
    auto output = out_info->Get(vtkDataObject::DATA_OBJECT());

    __check_not_null_ret(output, "No output object was found.");

    // Extract requested block
    if (this->BlockID < 0 || this->BlockID >= static_cast<int>(input->GetNumberOfBlocks()))
    {
        std::cerr << "Requested block does not exist." << std::endl;
        return 0;
    }

    auto block = input->GetBlock(this->BlockID);

    // Set output
    output->ShallowCopy(block);

    return 1;
}
