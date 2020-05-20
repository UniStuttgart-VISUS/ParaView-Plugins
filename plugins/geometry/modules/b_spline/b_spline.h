#pragma once

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyDataAlgorithm.h"

#include "Eigen/Dense"

#include <iterator>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

class b_spline : public vtkPolyDataAlgorithm
{
public:
    static b_spline *New();
    vtkTypeMacro(b_spline, vtkPolyDataAlgorithm);

    vtkGetMacro(NumberOfPoints, int);
    vtkSetMacro(NumberOfPoints, int);

    vtkGetMacro(Degree, int);
    vtkSetMacro(Degree, int);

    vtkGetMacro(HitEndpoints, int);
    vtkSetMacro(HitEndpoints, int);

protected:
    b_spline();

    int FillInputPortInformation(int, vtkInformation*);
    int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
    b_spline(const b_spline&);
    void operator=(const b_spline&);

    /// Compute point on the B-Spline
    Eigen::Vector3d compute_point(const std::vector<Eigen::Vector3d>& de_boor_points,
        std::vector<double>::const_iterator knot_vector_begin, std::vector<double>::const_iterator knot_vector_end,
        std::size_t degree, double arc_parameter) const;

    /// B-Spline basis function
    double basis_function(std::vector<double>::const_iterator knot_vector_begin, std::vector<double>::const_iterator knot_vector_end,
        std::size_t de_boor_index, std::size_t degree, double u) const;

    /// Create new control points representing the derivative of the B-Spline as new B-Spline
    std::vector<Eigen::Vector3d> derive(const std::vector<Eigen::Vector3d>& de_boor_points,
        std::vector<double>::const_iterator knot_vector_begin, std::vector<double>::const_iterator knot_vector_end, std::size_t degree) const;

    /// Number of output points
    int NumberOfPoints;

    /// Degree of the B-Spline
    int Degree;

    /// Hit endpoints through multiplicity in the knot vector
    int HitEndpoints;
};

/// Convert to string
template <typename t>
struct to_string
{
    std::string operator()(t value)
    {
        std::stringstream ss;
        ss << value;

        return ss.str();
    }
};

template <typename u>
struct to_string<std::pair<u, u>>
{
    std::string operator()(const std::pair<u, u>& value)
    {
        std::stringstream ss;
        ss << "[" << to_string<u>()(value.first) << ", " << to_string<u>()(value.second) << "]";

        return ss.str();
    }
};

/// Print a collection
template <typename it_t>
std::string print_collection(it_t begin, it_t end)
{
    // Concatenate all values
    std::stringstream ss;
    ss << "{";

    for (; begin != end; ++begin)
    {
        ss << " " << to_string<typename std::remove_reference<decltype(*begin)>::type>()(*begin) << ", ";
    }

    ss << "}";

    std::string string = ss.str();

    // Remove last comma
    const auto pos = string.find_last_of(',');

    if (pos != std::string::npos)
    {
        string.erase(pos, 1);
    }

    return string;
}

/// Get value at requested position from random access iterator
template <typename it_t>
typename std::iterator_traits<it_t>::reference access_at(it_t begin, std::size_t index)
{
    return *(begin + index);
}
