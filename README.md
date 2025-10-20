# Add Altitude to BRGM Geological Data

This project provides a Python script to enhance BRGM (Bureau de Recherches Géologiques et Minières) geological data by adding altitude information. The script supports multiple input formats: GeoJSON files, ESRI Shapefiles (.shp), and BRGM standard ZIP archives containing shapefiles. The altitude data is retrieved from the IGN (Institut Géographique National) France's elevation service.

## Features

- **Multiple Input Formats**: Supports GeoJSON, ESRI Shapefiles, and BRGM ZIP archives
- **Automatic Downloads**: Can download all French departments from BRGM website
- **Persistent Caching**: Saves downloaded data to avoid re-fetching during development
- **Format Preservation**: Output maintains the same format as input
- **Optimized API Usage**: 10km tiles with intelligent caching reduce network requests

## Prerequisites

- Python 3.12.7 (Anaconda distribution recommended)
- Required libraries: `requests`, `numpy`, `geopandas`
- Install with: `pip install requests numpy geopandas`

## Installation

1. Ensure Python and pip are installed and accessible via `python` and `pip` commands.
2. Install the required dependencies:
   ```
   pip install requests numpy geopandas
   ```

## Usage

The script supports two main modes: processing local files or downloading BRGM data.

### Process Local Files

```
python add_altitude.py --in <source_directory> --out <destination_directory> [--temp-dir <temp_directory>]
```

- `<source_directory>`: Directory containing input files (.geojson, .shp, .zip)
- `<destination_directory>`: Where processed files will be saved
- `<temp_directory>`: Optional persistent cache directory

### Download and Process BRGM Data

```
python add_altitude.py --from-brgm --out <destination_directory> [--departments <codes>] [--temp-dir <temp_directory>]
```

- `<destination_directory>`: Where downloaded and processed ZIPs will be saved
- `<codes>`: Comma-separated department codes (e.g., "073,074") or "all" for all departments
- `<temp_directory>`: Optional persistent cache directory

### Examples

```bash
# Process local files
python add_altitude.py --in ./brgm --out ./output --temp-dir ./temp

# Download and process all departments
python add_altitude.py --from-brgm --out ./output --departments all --temp-dir ./temp

# Download specific departments
python add_altitude.py --from-brgm --out ./output --departments 073,074 --temp-dir ./temp
```

## How It Works

1. **Input Processing**:
   - **GeoJSON files**: Load, process Point features, add `ALTITUDE` property, save as GeoJSON
   - **Shapefiles (.shp)**: Load with geopandas, filter Points, add `ALTITUDE` column, save as Shapefile
   - **BRGM ZIP files**: Extract shapefiles, process each .shp, re-compress with modified files

2. **Altitude Retrieval**:
   - For each point, calculate containing 10km × 10km tile
   - Check persistent cache first, download from IGN WMS if needed
   - Extract elevation value from tile pixel
   - Add to feature/geometry

3. **Caching System**:
   - **Memory cache**: Tiles loaded in RAM during session
   - **Disk cache**: Tiles saved as .npy files when `--temp-dir` specified
   - **ZIP cache**: Downloaded BRGM ZIPs kept when `--temp-dir` specified

### Supported Input Formats

- **GeoJSON**: FeatureCollections with Point geometries in EPSG:2154
- **ESRI Shapefile**: .shp files with Point geometries (other types ignored)
- **BRGM ZIP**: Standard format with shapefiles and documentation

### Department Codes

French departments use 3-digit codes (001-095), excluding 020 (replaced by 02A and 02B for Corsica).

## API Details

- **IGN WMS Service**: `https://data.geopf.fr/wms-r`
- **Layer**: `ELEVATION.ELEVATIONGRIDCOVERAGE.HIGHRES` (BDAlti 25m)
- **CRS**: `EPSG:2154` (Lambert 93)
- **Format**: `image/x-bil;bits=32`

Tiles are 10km × 10km (400 × 400 pixels at 25m resolution).

## Output

- **GeoJSON**: Additional `ALTITUDE` property (meters above sea level)
- **Shapefile**: Additional `ALTITUDE` column (float)
- **ZIP**: Re-compressed with modified shapefiles inside

Failed altitude queries set `ALTITUDE` to `null`/`NaN`.

## Error Handling

- Invalid directories or files logged and skipped
- API failures logged, processing continues
- Non-Point geometries ignored with warnings
- Network timeouts handled gracefully

## Performance Tips

- Use `--temp-dir` for development to avoid re-downloads
- Process related data together to maximize cache hits
- For large datasets, consider processing in batches

## License

This project is for educational purposes. Ensure compliance with BRGM and IGN data usage policies.