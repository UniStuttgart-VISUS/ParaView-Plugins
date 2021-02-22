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

    vtkSetStringMacro(RotationAngle);
    vtkGetStringMacro(RotationAngle);

    vtkSetStringMacro(RotationalVelocity);
    vtkGetStringMacro(RotationalVelocity);

    vtkDataArraySelection* GetScalarArraySelection();
    vtkDataArraySelection* GetVectorArraySelection();
    vtkDataArraySelection* GetVelocityArraySelection();
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

    char* RotationAngle;
    char* RotationalVelocity;

    vtkSmartPointer<vtkDataArraySelection> scalar_array_selection;
    vtkSmartPointer<vtkDataArraySelection> vector_array_selection;
    vtkSmartPointer<vtkDataArraySelection> velocity_array_selection;
    vtkSmartPointer<vtkDataArraySelection> pass_point_array_selection;
    vtkSmartPointer<vtkDataArraySelection> pass_cell_array_selection;
};
