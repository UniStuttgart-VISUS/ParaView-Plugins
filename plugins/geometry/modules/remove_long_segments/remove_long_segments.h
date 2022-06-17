#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyDataAlgorithm.h"

class remove_long_segments : public vtkPolyDataAlgorithm
{
public:
    static remove_long_segments *New();
    vtkTypeMacro(remove_long_segments, vtkPolyDataAlgorithm);

    vtkSetMacro(LengthFactor, double);
    vtkGetMacro(LengthFactor, double);

    vtkSetMacro(Absolute, int);
    vtkGetMacro(Absolute, int);

protected:
    remove_long_segments();

    int FillInputPortInformation(int, vtkInformation*);
    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    remove_long_segments(const remove_long_segments&);
    void operator=(const remove_long_segments&);

    /// Parameters
    double LengthFactor;
    int Absolute;
};
