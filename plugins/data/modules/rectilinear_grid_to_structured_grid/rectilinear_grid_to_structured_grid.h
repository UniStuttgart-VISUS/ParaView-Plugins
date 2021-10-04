#pragma once

#include "vtkStructuredGridAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT rectilinear_grid_to_structured_grid : public vtkStructuredGridAlgorithm
{
public:
    static rectilinear_grid_to_structured_grid* New();
    vtkTypeMacro(rectilinear_grid_to_structured_grid, vtkStructuredGridAlgorithm);

protected:
    rectilinear_grid_to_structured_grid();
    ~rectilinear_grid_to_structured_grid();

    int FillInputPortInformation(int, vtkInformation*);

    int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

    int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
    rectilinear_grid_to_structured_grid(const rectilinear_grid_to_structured_grid&);
    void operator=(const rectilinear_grid_to_structured_grid&);
};
