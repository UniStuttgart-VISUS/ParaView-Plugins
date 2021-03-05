#pragma once

#include "vtkRectilinearGridAlgorithm.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"

#include <array>

class VTK_EXPORT sample_points_to_grid : public vtkRectilinearGridAlgorithm
{
public:
    static sample_points_to_grid* New();
    vtkTypeMacro(sample_points_to_grid, vtkRectilinearGridAlgorithm);

    vtkSetVector3Macro(NumberOfBins, int);

protected:
    sample_points_to_grid();
    ~sample_points_to_grid();

    int FillInputPortInformation(int, vtkInformation*);

    int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

    int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
    sample_points_to_grid(const sample_points_to_grid&);
    void operator=(const sample_points_to_grid&);

    std::array<int, 3> NumberOfBins;
};
