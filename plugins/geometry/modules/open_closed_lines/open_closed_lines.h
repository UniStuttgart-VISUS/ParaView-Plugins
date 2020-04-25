#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyDataAlgorithm.h"

class open_closed_lines : public vtkPolyDataAlgorithm
{
public:
    static open_closed_lines *New();
    vtkTypeMacro(open_closed_lines, vtkPolyDataAlgorithm);

    vtkSetMacro(Method, int);
    vtkGetMacro(Method, int);

protected:
    open_closed_lines();

    int FillInputPortInformation(int, vtkInformation*);
    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    open_closed_lines(const open_closed_lines&);
    void operator=(const open_closed_lines&);

    /// Parameters
    int Method;
};
