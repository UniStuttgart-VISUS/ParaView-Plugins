# Average over Block

Average selected data array over all blocks from a multi block dataset, presenting it in its native data type.

## Input

The following inputs can be connected to the filter:

| Input      | Description                                           | Type               | Remark |
| ---------- | ----------------------------------------------------- | ------------------ | ------ |
| MultiBlock | Multiblock dataset with blocks over which to average. | Multiblock dataset |        |

## Parameters

The following parameters are available in the properties panel in ParaView:

| Parameter  | Description                                                  | Default value |
| ---------- | ------------------------------------------------------------ | ------------- |
| Data array | Name of the data array used for averaging (component-wise). The number of elements in each block must match. |               |

## Output

The output is the geometry of the first input block, attached with the array containing the averaged data values.
