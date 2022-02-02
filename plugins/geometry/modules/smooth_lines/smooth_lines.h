#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSetAlgorithm.h"
#include "vtkStructuredGrid.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"

class smooth_lines : public vtkMultiBlockDataSetAlgorithm
{
public:
    static smooth_lines *New();
    vtkTypeMacro(smooth_lines, vtkMultiBlockDataSetAlgorithm);

    vtkSetMacro(NumIterations, int);
    vtkGetMacro(NumIterations, int);

    vtkSetMacro(Method, int);
    vtkGetMacro(Method, int);

    vtkSetMacro(Variant, int);
    vtkGetMacro(Variant, int);

    vtkSetMacro(ImplicitLambda, double);
    vtkGetMacro(ImplicitLambda, double);

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

    /// Enums
    enum class method_t
    {
        implicit_gaussian, gaussian, taubin, perp_gaussian
    };

    enum class variant_t
    {
        normal, fixed_endpoints, growing
    };

    /*
     * Implicit Gaussian smoothing
     *
     * @param vertices Line points
     * @param matrix Matrix for solving the implicit problem
     *
     * @return Displaced points
     */
    void gaussian_smoothing(Eigen::MatrixXd& vertices, const Eigen::SparseMatrix<double>& matrix) const;

    /*
     * Gaussian smoothing
     *
     * @param points Line points
     * @param index Index of the point for which the displacement from smoothing is calculated
     * @param weight Smoothing weight; positive for smoothing, negative for inflation
     * @param variant Variant of the smoothing method
     *
     * @return Displaced point
     */
    Eigen::Vector3d gaussian_smoothing(const std::vector<Eigen::Vector3d>& points,
        std::size_t index, double weight, variant_t variant) const;

    /*
     * Gaussian smoothing restricted to movement along the given vector field
     *
     * @param points Line points
     * @param index Index of the point for which the displacement from smoothing is calculated
     * @param weight Smoothing weight; positive for smoothing, negative for inflation
     *
     * @return Displaced point
     */
    Eigen::Vector3d gaussian_smoothing(const std::vector<Eigen::Vector3d>& points,
        std::size_t index, double weight, const vtkStructuredGrid* vector_field) const;

    /// Parameters
    int NumIterations;
    int Method;
    int Variant;
    double ImplicitLambda, Lambda, Mu;
};
