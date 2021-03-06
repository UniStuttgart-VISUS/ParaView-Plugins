#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyDataAlgorithm.h"

#include "Eigen/Dense"

class smooth_lines : public vtkPolyDataAlgorithm
{
public:
    static smooth_lines *New();
    vtkTypeMacro(smooth_lines, vtkPolyDataAlgorithm);

    vtkSetMacro(NumIterations, int);
    vtkGetMacro(NumIterations, int);

    vtkSetMacro(Lambda, double);
    vtkGetMacro(Lambda, double);

    vtkSetMacro(Mu, double);
    vtkGetMacro(Mu, double);

protected:
    smooth_lines();

    int FillInputPortInformation(int, vtkInformation*);
    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    smooth_lines(const smooth_lines&);
    void operator=(const smooth_lines&);

    /*
     * Gaussian smoothing
     *
     * @param points Line points
     * @param index Index of the point for which the displacement from smoothing is calculated
     * @param weight Smoothing weight; positive for smoothing, negative for inflation
     *
     * @return Displaced point
     */
    Eigen::Vector3d gaussian_smoothing(const std::vector<Eigen::Vector3d>& points, std::size_t index, double weight) const;

    /// Parameters
    int NumIterations;
    double Lambda, Mu;
};
