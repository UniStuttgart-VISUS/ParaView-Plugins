#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyDataAlgorithm.h"

class truncate_lines : public vtkPolyDataAlgorithm
{
public:
    static truncate_lines *New();
    vtkTypeMacro(truncate_lines, vtkPolyDataAlgorithm);

    vtkSetMacro(Offset, int);
    vtkGetMacro(Offset, int);

    vtkSetMacro(NumPoints, int);
    vtkGetMacro(NumPoints, int);

protected:
    truncate_lines();

    int FillInputPortInformation(int, vtkInformation*);
    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    truncate_lines(const truncate_lines&);
    void operator=(const truncate_lines&);

    /// Parameters
    int Offset;
    int NumPoints;
};
