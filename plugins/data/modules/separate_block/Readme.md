# Separate Block

Extract a block from a multi block dataset, presenting it in its native data type.

*Please note that the output datatype defaults to an unstructured grid, and therefore does not work for other datatypes when loading from a pvsm state.*

## Input

The following inputs can be connected to the filter:

| Input                     | Description                                                                               | Type                  | Remark        |
|---------------------------|-------------------------------------------------------------------------------------------|-----------------------|---------------|
| MultiBlock                | Multiblock dataset, from which to extract a block.                                        | Multiblock dataset    |               |

## Parameters

The following parameters are available in the properties panel in ParaView:

| Parameter                 | Description                                                                                                   | Default value         |
|---------------------------|---------------------------------------------------------------------------------------------------------------|-----------------------|
| Block ID                  | ID of the block to extract.                                                                                   | 0                     |
