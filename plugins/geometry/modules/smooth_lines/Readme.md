# Smooth lines

Smooth lines by applying Gaussian or Taubin smoothing. Additionally, a restricted smoothing can be applied, which only updates the positions perpendicular to a given vector field.

## Input

The following inputs can be connected to the filter:

| Input                     | Description                                                               | Type              | Remark        |
|---------------------------|---------------------------------------------------------------------------|-------------------|---------------|
| Lines                     | Lines which are going to be smoothed.                                     | Poly data         |               |
| Grid                      | Grid containing a vector field, used to restrict smoothing.               | Structured Grid   |               |

## Parameters

The following parameters are available in the properties panel in ParaView:

| Parameter                                                                             | Description                                                                           | Default value |
|---------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------|---------------|
| Vectors                                                                               | Vector field for restricted, perpendicular smoothing.                                 | 100           |
| MaxNumIterations                                                                      | Number of smoothing steps performed.                                                  | 100           |
| Method                                                                                | Different smoothing methods: Gaussian, Taubin, and perpendicular Gaussian.            | Gaussian      |
| Lambda (![Equation](https://render.githubusercontent.com/render/math?math=\lambda))   | Gaussian smoothing factor.                                                            | 0.7           |
| Mu (![Equation](https://render.githubusercontent.com/render/math?math=\mu))           | Inflation factor for Taubin smoothing (set zero for Gaussian smoothing).              | -0.75         |

Both smoothing parameters are values between 0 (exlusive; no smoothing) and 1 (inclusive; large smoothing). Note that for the inflation factor of Taubin smoothing, ![Equation](https://render.githubusercontent.com/render/math?math=\mu) has to meet the following additional constraint:

(1) ![Equation](https://render.githubusercontent.com/render/math?math=|\lambda|\\,\\,\leq\\,\\,-|\mu|).
