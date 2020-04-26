# Truncate lines

Truncate lines by defining an offset and a number of points for the new lines.

## Input

The following inputs can be connected to the filter:

| Input                     | Description                                                               | Type          | Remark        |
|---------------------------|---------------------------------------------------------------------------|---------------|---------------|
| Lines                     | Lines which are going to be truncated.                                    | Poly data     |               |

## Parameters

The following parameters are available in the properties panel in ParaView:

| Parameter         | Description                                                                           | Default value         |
|-------------------|---------------------------------------------------------------------------------------|-----------------------|
| Offset            | Offset defining the start index of the new lines.                                     | 0                     |
| Number of points  | Maximum number of points of the new line, where everything behind is truncated.       | 1000                  |
