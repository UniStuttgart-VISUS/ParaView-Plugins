#pragma once

#include "vtkRectilinearGridAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT structured_grid_to_rectilinear_grid : public vtkRectilinearGridAlgorithm
{
public:
    static structured_grid_to_rectilinear_grid* New();
    vtkTypeMacro(structured_grid_to_rectilinear_grid, vtkRectilinearGridAlgorithm);

protected:
    structured_grid_to_rectilinear_grid();
    ~structured_grid_to_rectilinear_grid();

    int FillInputPortInformation(int, vtkInformation*);

    int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

    int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
    structured_grid_to_rectilinear_grid(const structured_grid_to_rectilinear_grid&);
    void operator=(const structured_grid_to_rectilinear_grid&);
};
