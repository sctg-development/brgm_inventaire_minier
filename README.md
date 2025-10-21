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

## Automated Release Process

All French departments are automatically processed and released using GitHub Actions workflows. The `release.yml` workflow:

1. **Parallel Processing**: Processes all 13 metropolitan French regions in parallel
2. **Regional Organization**: Groups departments logically by administrative region
3. **Automatic Releases**: Creates weekly releases with all processed data
4. **Individual Downloads**: Each region's data is available as a separate ZIP file

### How to Download Files

1. **Go to Releases**: Visit the project's [Releases page](https://github.com/sctg-development/brgm_inventaire_minier/releases)
2. **Select Region**: Each release contains ZIP files for all 13 regions:
   - `GEO050K_HARM_*.zip` - One file per department (except 59-62 combined, and 75-77-78-91-92-93-94-95 combined)
3. **Download ZIP**: Click the department's ZIP file to download
4. **Verify Integrity**: Each ZIP includes a build provenance attestation

### Verifying Attestations

GitHub provides SLSA provenance attestations for each file to verify integrity and origin:

1. **Download Attestation**: Each ZIP file's attestation is generated during release
2. **Verify with GitHub CLI**:
   ```bash
   gh attestation verify GEO050K_HARM_075_077_078_091_092_093_094_095.zip --owner sctg-development
   ```
3. **Requirements**: GitHub CLI (`gh`) must be installed and authenticated

### Coverage by Release

Each release contains data for all 96 French departments across 13 regions:

- **Auvergne-Rhône-Alpes**: 12 departments (001, 003, 007, 015, 026, 038, 042, 043, 063, 069, 073, 074)
- **Bourgogne-Franche-Comté**: 8 departments (021, 025, 039, 058, 070, 071, 089, 090)
- **Bretagne**: 4 departments (022, 029, 035, 056)
- **Centre-Val de Loire**: 6 departments (018, 028, 036, 037, 041, 045)
- **Corse**: 2 departments (02A, 02B)
- **Grand Est**: 10 departments (008, 010, 051, 052, 054, 055, 057, 067, 068, 088)
- **Hauts-de-France**: 5 departments (002, 059, 060, 062, 080)
- **Île-de-France**: 8 departments (075, 077, 078, 091, 092, 093, 094, 095)
- **Normandie**: 5 departments (014, 027, 050, 061, 076)
- **Nouvelle-Aquitaine**: 12 departments (016, 017, 019, 023, 024, 033, 040, 047, 064, 079, 086, 087)
- **Occitanie**: 12 departments (009, 011, 030, 031, 032, 034, 046, 048, 065, 066, 081, 082)
- **Pays de la Loire**: 5 departments (044, 049, 053, 072, 085)
- **Provence-Alpes-Côte d'Azur**: 6 departments (004, 005, 006, 013, 083, 084)

## License

This project is licensed under the MIT License for educational purposes. Ensure compliance with BRGM and IGN data usage policies.