#pragma once

#include "vtkImageAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

class VTK_EXPORT rectilinear_grid_to_image : public vtkImageAlgorithm
{
public:
    static rectilinear_grid_to_image* New();
    vtkTypeMacro(rectilinear_grid_to_image, vtkImageAlgorithm);

protected:
    rectilinear_grid_to_image();
    ~rectilinear_grid_to_image();

    int FillInputPortInformation(int, vtkInformation*);

    int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

    int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
    rectilinear_grid_to_image(const rectilinear_grid_to_image&);
    void operator=(const rectilinear_grid_to_image&);
};
