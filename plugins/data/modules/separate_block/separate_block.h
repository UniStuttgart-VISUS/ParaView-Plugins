#pragma once

#include "vtkAlgorithm.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"

class separate_block : public vtkAlgorithm
{
public:
    static separate_block *New();
    vtkTypeMacro(separate_block, vtkAlgorithm);

    vtkGetMacro(BlockID, int);
    vtkSetMacro(BlockID, int);

protected:
    separate_block();

    virtual int FillInputPortInformation(int, vtkInformation*) override;
    virtual int FillOutputPortInformation(int, vtkInformation*) override;

    virtual int ProcessRequest(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

    virtual int RequestDataObject(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
    virtual int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
    virtual int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
    virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
    separate_block(const separate_block&);
    void operator=(const separate_block&);

    /// ID of the block to extract
    int BlockID;
};
