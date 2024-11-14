# Source



# Download Steps

## Auto download

The package provides a function: `SeisTools.Topography.download_srtm_block(lat::Real, lon::Real)`. You can download the tile covering the specified `lat` and `lon`. The `lat` and `lon` do not need to be the center of tile.

## Manual download

1. Open the [link to SRTM90 download page](https://srtm.csi.cgiar.org/srtmdata/)
   ![homepage of SRTM90](figs/page1.png)
2. choose Tile 5x5 degree block, `Esri ASCII` format, then click on blocks in map to select
3. click search to list results
   ![search result](figs/page2.png)
4. click `Download SRTM` to download data of this block and save it into `.../SeisTools.jl/external/Topography/SRTM90/raw_data`
