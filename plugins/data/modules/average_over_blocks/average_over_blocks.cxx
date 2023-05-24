#include "average_over_blocks.h"

#include "common/checks.h"

#include "vtkDataObject.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkMultiPieceDataSet.h"
#include "vtkObjectFactory.h"
#include "vtkCellData.h"
#include "vtkFieldData.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkRectilinearGrid.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStructuredGrid.h"
#include "vtkStructuredPoints.h"
#include "vtkType.h"
#include "vtkUnstructuredGrid.h"

#include <functional>
#include <iostream>
#include <string>

vtkStandardNewMacro(average_over_blocks);

average_over_blocks::average_over_blocks()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

int average_over_blocks::ProcessRequest(vtkInformation* request, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
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

int average_over_blocks::RequestDataObject(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get input block
    auto in_info = input_vector[0]->GetInformationObject(0);
    auto input = vtkMultiBlockDataSet::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(input, "Input is not a multiblock dataset.");

    // If there is no block or the specified block does not exist, create a dummy object
    // This step is necessary to allow loading from a state
    if (input->GetNumberOfBlocks() == 0)
    {
        return 0;
    }

    // Create appropriate data object if the type differs from the previous execution
    auto output = output_vector->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT());

    if (output == nullptr || output->GetDataObjectType() != input->GetBlock(0)->GetDataObjectType())
    {
        // Create output object of appropriate type
        auto* newOutput = input->GetBlock(0)->NewInstance();
        output_vector->GetInformationObject(0)->Set(vtkDataObject::DATA_OBJECT(), newOutput);
        newOutput->Delete();
    }

    return 1;
}

int average_over_blocks::RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector* output_vector)
{
    return 1;
}

int average_over_blocks::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
        return 1;
    }

    return 0;
}

int average_over_blocks::FillOutputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataObject");
        return 1;
    }

    return 0;
}

int average_over_blocks::RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*)
{
    return 1;
}

int average_over_blocks::RequestData(vtkInformation*, vtkInformationVector** input_vector, vtkInformationVector* output_vector)
{
    // Get access to information and data
    auto in_info = input_vector[0]->GetInformationObject(0);
    auto input = vtkMultiBlockDataSet::SafeDownCast(in_info->Get(vtkDataObject::DATA_OBJECT()));

    __check_not_null_ret(input, "Input is not a multiblock dataset.");

    // Get output
    auto out_info = output_vector->GetInformationObject(0);
    auto output = out_info->Get(vtkDataObject::DATA_OBJECT());

    __check_not_null_ret(output, "No output object was found.");

    // Copy first block
    if (input->GetNumberOfBlocks() == 0)
    {
        return 1;
    }

    output->DeepCopy(input->GetBlock(0));

    // Handle multi-piece datasets
    std::function<vtkDataSet* (vtkMultiBlockDataSet*, int, int)> get_input_data;
    std::function<vtkDataSet* (vtkDataObject*, int)> get_output_data;

    vtkIdType num_pieces = 1;

    if (vtkMultiPieceDataSet::SafeDownCast(input->GetBlock(0)))
    {
        get_input_data = [](vtkMultiBlockDataSet* input, int block, int piece)
        {
            return vtkMultiPieceDataSet::SafeDownCast(input->GetBlock(block))->GetPiece(piece);
        };

        get_output_data = [](vtkDataObject* output, int piece)
        {
            return vtkMultiPieceDataSet::SafeDownCast(output)->GetPiece(piece);
        };

        num_pieces = vtkMultiPieceDataSet::SafeDownCast(input->GetBlock(0))->GetNumberOfPieces();
    }
    else
    {
        get_input_data = [](vtkMultiBlockDataSet* input, int block, int piece)
        {
            return vtkDataSet::SafeDownCast(input->GetBlock(block));
        };

        get_output_data = [](vtkDataObject* output, int piece)
        {
            return vtkDataSet::SafeDownCast(output);
        };
    }

    // Average array over all blocks
    auto pdGetter = [](vtkDataSet* ds, const std::string& name) { return ds->GetPointData()->GetArray(name.c_str()); };
    auto cdGetter = [](vtkDataSet* ds, const std::string& name) { return ds->GetCellData()->GetArray(name.c_str()); };
    auto fdGetter = [](vtkDataSet* ds, const std::string& name) { return ds->GetFieldData()->GetArray(name.c_str()); };

    std::function<vtkDataArray* (vtkDataSet*)> getter;

    int association = vtkDataObject::FIELD_ASSOCIATION_NONE;
    const std::string name = GetInputArrayToProcess(0, get_input_data(input, 0, 0), association)->GetName();

    if (association == vtkDataObject::FIELD_ASSOCIATION_POINTS)
    {
        getter = std::bind(pdGetter, std::placeholders::_1, name);
    }
    else if (association == vtkDataObject::FIELD_ASSOCIATION_CELLS)
    {
        getter = std::bind(cdGetter, std::placeholders::_1, name);
    }
    else
    {
        getter = std::bind(fdGetter, std::placeholders::_1, name);
    }

    for (vtkIdType p = 0; p < num_pieces; ++p)
    {
        auto* average = getter(get_output_data(output, p));

        if (!average) {
            std::cerr << "Data array must exist for all blocks (missing for: 0)." << std::endl;
            return 0;
        }

        for (vtkIdType b = 1; b < input->GetNumberOfBlocks(); ++b)
        {
            auto* data = getter(get_input_data(input, b, p));

            if (!data) {
                std::cerr << "Data array must exist for all blocks (missing for: " << b << ")." << std::endl;
                return 0;
            }

            if (data->GetNumberOfValues() != average->GetNumberOfValues())
            {
                std::cerr << "Number of elements does not match for all blocks." << std::endl;
                return 0;
            }

            for (vtkIdType i = 0; i < data->GetNumberOfValues(); ++i)
            {
                average->SetVariantValue(i, average->GetVariantValue(i).ToDouble() + data->GetVariantValue(i).ToDouble());
            }
        }

        for (vtkIdType i = 0; i < average->GetNumberOfValues(); ++i)
        {
            average->SetVariantValue(i, average->GetVariantValue(i).ToDouble() / input->GetNumberOfBlocks());
        }
    }

    return 1;
}
