#pragma once

#include "vtkAlgorithm.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"

#include <array>

class VTK_EXPORT regular_grid : public vtkAlgorithm
{
    public:
        static regular_grid* New();
        vtkTypeMacro(regular_grid, vtkAlgorithm);

        vtkSetMacro(OutputType, int);
        vtkSetMacro(Anchor, int);

        vtkSetVector3Macro(Dimension, int);
        vtkSetVector3Macro(Origin, double);
        vtkSetVector3Macro(CellSize, double);

    protected:
        regular_grid();

        virtual int FillOutputPortInformation(int, vtkInformation*) override;

        virtual int ProcessRequest(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

        virtual int RequestDataObject(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
        virtual int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
        virtual int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
        virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

    private:
        regular_grid(const regular_grid&);
        void operator=(const regular_grid&);

        /// Output type
        int OutputType;

        /// Anchor for the origin
        int Anchor;

        /// Extent of the grid in number of cells
        std::array<int, 3> Dimension;

        /// Origin
        std::array<double, 3> Origin;

        /// Cell size
        std::array<double, 3> CellSize;
};
