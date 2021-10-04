#pragma once

#include "vtkUnstructuredGridAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT structured_grid_to_unstructured_grid : public vtkUnstructuredGridAlgorithm
{
public:
    static structured_grid_to_unstructured_grid* New();
    vtkTypeMacro(structured_grid_to_unstructured_grid, vtkUnstructuredGridAlgorithm);

protected:
    structured_grid_to_unstructured_grid();
    ~structured_grid_to_unstructured_grid();

    int FillInputPortInformation(int, vtkInformation*);

    int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

    int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
    structured_grid_to_unstructured_grid(const structured_grid_to_unstructured_grid&);
    void operator=(const structured_grid_to_unstructured_grid&);
};
