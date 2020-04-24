#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyDataAlgorithm.h"

class connect_points : public vtkPolyDataAlgorithm
{
public:
    static connect_points *New();
    vtkTypeMacro(connect_points, vtkPolyDataAlgorithm);

protected:
    connect_points();

    int FillInputPortInformation(int, vtkInformation*);
    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    connect_points(const connect_points&);
    void operator=(const connect_points&);
};
