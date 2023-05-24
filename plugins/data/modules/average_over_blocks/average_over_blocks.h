#pragma once

#include "vtkAlgorithm.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"

class average_over_blocks : public vtkAlgorithm
{
public:
    static average_over_blocks *New();
    vtkTypeMacro(average_over_blocks, vtkAlgorithm);

protected:
    average_over_blocks();

    virtual int FillInputPortInformation(int, vtkInformation*) override;
    virtual int FillOutputPortInformation(int, vtkInformation*) override;

    virtual int ProcessRequest(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

    virtual int RequestDataObject(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
    virtual int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
    virtual int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
    virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
    average_over_blocks(const average_over_blocks&);
    void operator=(const average_over_blocks&);
};
