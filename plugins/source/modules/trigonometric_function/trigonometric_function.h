#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyDataAlgorithm.h"

class trigonometric_function : public vtkPolyDataAlgorithm
{
public:
    static trigonometric_function *New();
    vtkTypeMacro(trigonometric_function, vtkPolyDataAlgorithm);

    vtkSetMacro(NumPoints, int);
    vtkGetMacro(NumPoints, int);

    vtkSetMacro(Length, double);
    vtkGetMacro(Length, double);

    vtkSetMacro(Function, int);
    vtkGetMacro(Function, int);

    vtkSetMacro(Alpha, double);
    vtkGetMacro(Alpha, double);

    vtkSetMacro(Beta, double);
    vtkGetMacro(Beta, double);

    vtkSetMacro(Gamma, double);
    vtkGetMacro(Gamma, double);

    vtkSetMacro(Delta, double);
    vtkGetMacro(Delta, double);

protected:
    trigonometric_function();

    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    trigonometric_function(const trigonometric_function&);
    void operator=(const trigonometric_function&);

    int NumPoints;
    double Length;

    int Function;

    double Alpha;
    double Beta;
    double Gamma;
    double Delta;
};
