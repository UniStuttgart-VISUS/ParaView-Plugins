#pragma once

#include "vtkRectilinearGridAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT fix_errors_in_grids : public vtkRectilinearGridAlgorithm
{
public:
    static fix_errors_in_grids* New();
    vtkTypeMacro(fix_errors_in_grids, vtkRectilinearGridAlgorithm);

protected:
    fix_errors_in_grids();
    ~fix_errors_in_grids();

    int FillInputPortInformation(int, vtkInformation*);

    int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

    int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
    fix_errors_in_grids(const fix_errors_in_grids&);
    void operator=(const fix_errors_in_grids&);
};
