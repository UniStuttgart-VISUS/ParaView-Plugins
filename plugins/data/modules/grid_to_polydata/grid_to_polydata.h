#pragma once

#include "vtkPolyDataAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT grid_to_polydata : public vtkPolyDataAlgorithm
{
public:
    static grid_to_polydata* New();
    vtkTypeMacro(grid_to_polydata, vtkPolyDataAlgorithm);

protected:
    grid_to_polydata();
    ~grid_to_polydata();

    int FillInputPortInformation(int, vtkInformation*);

    int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

    int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
    grid_to_polydata(const grid_to_polydata&);
    void operator=(const grid_to_polydata&);
};
