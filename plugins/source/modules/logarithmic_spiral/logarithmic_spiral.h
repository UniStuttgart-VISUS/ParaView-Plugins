#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyDataAlgorithm.h"

class logarithmic_spiral : public vtkPolyDataAlgorithm
{
public:
    static logarithmic_spiral *New();
    vtkTypeMacro(logarithmic_spiral, vtkPolyDataAlgorithm);

    vtkSetMacro(NumPoints, int);
    vtkGetMacro(NumPoints, int);

    vtkSetMacro(Length, double);
    vtkGetMacro(Length, double);

    vtkSetMacro(SizeFactor, double);
    vtkGetMacro(SizeFactor, double);

    vtkSetMacro(Slope, double);
    vtkGetMacro(Slope, double);

protected:
    logarithmic_spiral();

    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    logarithmic_spiral(const logarithmic_spiral&);
    void operator=(const logarithmic_spiral&);

    int NumPoints;
    double Length;

    double SizeFactor;
    double Slope;
};
