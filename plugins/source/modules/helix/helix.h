#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyDataAlgorithm.h"

class helix : public vtkPolyDataAlgorithm
{
public:
    static helix *New();
    vtkTypeMacro(helix, vtkPolyDataAlgorithm);

    vtkSetMacro(NumPoints, int);
    vtkGetMacro(NumPoints, int);

    vtkSetMacro(Length, double);
    vtkGetMacro(Length, double);

    vtkSetMacro(Radius, double);
    vtkGetMacro(Radius, double);

    vtkSetMacro(Windings, double);
    vtkGetMacro(Windings, double);

protected:
    helix();

    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    helix(const helix&);
    void operator=(const helix&);

    double NumPoints;

    double Length;
    double Radius;
    double Windings;
};
