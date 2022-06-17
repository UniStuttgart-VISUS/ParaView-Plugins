#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyDataAlgorithm.h"

class dense_grid : public vtkPolyDataAlgorithm
{
public:
    static dense_grid *New();
    vtkTypeMacro(dense_grid, vtkPolyDataAlgorithm);

    vtkSetMacro(Ratio, int);
    vtkGetMacro(Ratio, int);

protected:
    dense_grid();

    int FillInputPortInformation(int, vtkInformation*);
    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    dense_grid(const dense_grid&);
    void operator=(const dense_grid&);

    /// Parameters
    int Ratio;
};
