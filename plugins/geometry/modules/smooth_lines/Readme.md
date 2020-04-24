# Smooth lines

Smooth lines by applying Gaussian or Taubin smoothing.

## Input

The following inputs can be connected to the filter:

| Input                     | Description                                                               | Type          | Remark        |
|---------------------------|---------------------------------------------------------------------------|---------------|---------------|
| Lines                     | Lines which are going to be smoothed.                                     | Poly data     |               |

## Parameters

The following parameters are available in the properties panel in ParaView:

| Parameter                                                                             | Description                                                                           | Default value |
|---------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------|---------------|
| MaxNumIterations                                                                      | Number of smoothing steps performed.                                                  | 100           |
| Lambda (![Equation](https://render.githubusercontent.com/render/math?math=\lambda))   | Gaussian smoothing factor.                                                            | 0.7           |
| Mu (![Equation](https://render.githubusercontent.com/render/math?math=\mu))           | Inflation factor for Taubin smoothing (set zero for Gaussian smoothing).              | -0.75         |

Both smoothing parameters are values between 0 (no smoothing) and 1 (large smoothing). Note that for non-zero values for the inflaction factor, ![Equation](https://render.githubusercontent.com/render/math?math=\mu) has to meet the following additional constraint:

(1) ![Equation](https://render.githubusercontent.com/render/math?math=|\lambda|\\,\\,\leq\\,\\,-|\mu|).
