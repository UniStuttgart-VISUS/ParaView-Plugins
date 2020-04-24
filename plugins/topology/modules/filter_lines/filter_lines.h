#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyDataAlgorithm.h"

class filter_lines : public vtkPolyDataAlgorithm
{
public:
    static filter_lines *New();
    vtkTypeMacro(filter_lines, vtkPolyDataAlgorithm);

    vtkSetMacro(AngleFilter, int);
    vtkGetMacro(AngleFilter, int);

    vtkSetMacro(MaxAngle, double);
    vtkGetMacro(MaxAngle, double);

    vtkSetMacro(SizeFilter, int);
    vtkGetMacro(SizeFilter, int);

    vtkSetMacro(MinSize, int);
    vtkGetMacro(MinSize, int);

protected:
    filter_lines();

    int FillInputPortInformation(int, vtkInformation*);
    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    filter_lines(const filter_lines&);
    void operator=(const filter_lines&);

    /// Parameters
    int AngleFilter;
    double MaxAngle;

    int SizeFilter;
    int MinSize;
};
