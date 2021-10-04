#pragma once

#include "vtkRectilinearGridAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT image_to_rectilinear_grid : public vtkRectilinearGridAlgorithm
{
public:
    static image_to_rectilinear_grid* New();
    vtkTypeMacro(image_to_rectilinear_grid, vtkRectilinearGridAlgorithm);

protected:
    image_to_rectilinear_grid();
    ~image_to_rectilinear_grid();

    int FillInputPortInformation(int, vtkInformation*);

    int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

    int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
    image_to_rectilinear_grid(const image_to_rectilinear_grid&);
    void operator=(const image_to_rectilinear_grid&);
};
