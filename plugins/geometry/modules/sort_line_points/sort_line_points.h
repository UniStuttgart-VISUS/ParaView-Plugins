#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyDataAlgorithm.h"

class sort_line_points : public vtkPolyDataAlgorithm
{
public:
    static sort_line_points *New();
    vtkTypeMacro(sort_line_points, vtkPolyDataAlgorithm);

protected:
    sort_line_points();

    int FillInputPortInformation(int, vtkInformation*);
    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    sort_line_points(const sort_line_points&);
    void operator=(const sort_line_points&);
};
