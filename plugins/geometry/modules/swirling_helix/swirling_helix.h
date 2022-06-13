#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyDataAlgorithm.h"

class swirling_helix : public vtkPolyDataAlgorithm
{
public:
    static swirling_helix *New();
    vtkTypeMacro(swirling_helix, vtkPolyDataAlgorithm);

    vtkSetMacro(Radius, double);
    vtkGetMacro(Radius, double);

    vtkSetMacro(Windings, double);
    vtkGetMacro(Windings, double);

protected:
    swirling_helix();

    int FillInputPortInformation(int, vtkInformation*);
    int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    swirling_helix(const swirling_helix&);
    void operator=(const swirling_helix&);

    double Radius;
    double Windings;
};
