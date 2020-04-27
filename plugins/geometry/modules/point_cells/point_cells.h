#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyDataAlgorithm.h"

class point_cells : public vtkPolyDataAlgorithm
{
public:
    static point_cells *New();
    vtkTypeMacro(point_cells, vtkPolyDataAlgorithm);

protected:
    point_cells();

    int FillInputPortInformation(int, vtkInformation*);
    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    point_cells(const point_cells&);
    void operator=(const point_cells&);
};
