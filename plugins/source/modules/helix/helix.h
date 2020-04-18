#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyDataAlgorithm.h"

class Helix : public vtkPolyDataAlgorithm
{
public:
    static Helix *New();
    vtkTypeMacro(Helix, vtkPolyDataAlgorithm);

    vtkSetMacro(NumPoints, int);
    vtkGetMacro(NumPoints, int);

    vtkSetMacro(Length, double);
    vtkGetMacro(Length, double);

    vtkSetMacro(Radius, double);
    vtkGetMacro(Radius, double);

    vtkSetMacro(Windings, double);
    vtkGetMacro(Windings, double);

protected:
    Helix();

    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    Helix(const Helix&);
    void operator=(const Helix&);

    double NumPoints;

    double Length;
    double Radius;
    double Windings;
};
