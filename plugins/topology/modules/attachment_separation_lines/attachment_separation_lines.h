#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyDataAlgorithm.h"

class attachment_separation_lines : public vtkPolyDataAlgorithm
{
public:
    static attachment_separation_lines *New();
    vtkTypeMacro(attachment_separation_lines, vtkPolyDataAlgorithm);

protected:
    attachment_separation_lines();

    int FillInputPortInformation(int, vtkInformation*);
    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    attachment_separation_lines(const attachment_separation_lines&);
    void operator=(const attachment_separation_lines&);
};
