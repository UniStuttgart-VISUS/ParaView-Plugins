#pragma once

#include "vtkRectilinearGridAlgorithm.h"

#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkSmartPointer.h"

class VTK_EXPORT resample_rotating_grid : public vtkRectilinearGridAlgorithm
{
public:
    static resample_rotating_grid* New();
    vtkTypeMacro(resample_rotating_grid, vtkRectilinearGridAlgorithm);

    vtkSetStringMacro(RotationColumn);
    vtkGetStringMacro(RotationColumn);

    vtkDataArraySelection* GetScalarArraySelection();
    vtkDataArraySelection* GetVectorArraySelection();
    vtkDataArraySelection* GetPassPointArraySelection();
    vtkDataArraySelection* GetPassCellArraySelection();

protected:
    resample_rotating_grid();
    ~resample_rotating_grid();

    int FillInputPortInformation(int, vtkInformation*);

    int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

    int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
    resample_rotating_grid(const resample_rotating_grid&);
    void operator=(const resample_rotating_grid&);

    char* RotationColumn;

    vtkSmartPointer<vtkDataArraySelection> scalar_array_selection;
    vtkSmartPointer<vtkDataArraySelection> vector_array_selection;
    vtkSmartPointer<vtkDataArraySelection> pass_point_array_selection;
    vtkSmartPointer<vtkDataArraySelection> pass_cell_array_selection;
};
