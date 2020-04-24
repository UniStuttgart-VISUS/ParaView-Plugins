#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyDataAlgorithm.h"

class connect_lines : public vtkPolyDataAlgorithm
{
public:
    static connect_lines *New();
    vtkTypeMacro(connect_lines, vtkPolyDataAlgorithm);

protected:
    connect_lines();

    int FillInputPortInformation(int, vtkInformation*);
    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    connect_lines(const connect_lines&);
    void operator=(const connect_lines&);
};
