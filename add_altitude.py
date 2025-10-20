#!/usr/bin/env python3
# Copyright (c) 2025 Ronan Le Meillat - SCTG Development

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""
add_altitude.py

This script adds altitude data to GeoJSON point features from BRGM geological data.
It reads GeoJSON files from a source directory, queries the IGN WMTS elevation service for each point,
adds the altitude to the feature properties, and saves the modified files to a destination directory.
(c) 2025 Ronan Le Meillat - SCTG Development

Usage:
    python add_altitude.py <source_directory> <destination_directory>

Arguments:
    source_directory: Path to the directory containing GeoJSON files.
    destination_directory: Path to the directory where modified GeoJSON files will be saved.

Dependencies:
    - requests: For making HTTP requests.
    - owslib: For accessing WMTS services.
    - pyproj: For coordinate reprojection.
    - rasterio: For reading raster data from tiles.
    - Install with: pip install requests owslib pyproj rasterio

Note: This script assumes all GeoJSON files contain point geometries in EPSG:2154 (Lambert 93).
"""

import os
import sys
import json
import requests
import numpy as np
import time
import zipfile
import tempfile
import shutil
from pathlib import Path
import geopandas as gpd

# IGN WMS service URL
WMS_URL = 'https://data.geopf.fr/wms-r'
LAYER = 'ELEVATION.ELEVATIONGRIDCOVERAGE.HIGHRES'

# Tile parameters for caching
TILE_SIZE = 10000  # 10km tiles
PIXEL_SIZE = 25    # 25m resolution
TILE_PIXELS = TILE_SIZE // PIXEL_SIZE  # 400 pixels

# Cache for downloaded tiles
tile_cache = {}

# Temp dir for persistent cache
persistent_temp_dir = None

def download_tile_with_retry(url, max_retries=3, delay=60):
    """
    Download a tile with retry logic for network errors.
    
    Args:
        url (str): The URL to download from
        max_retries (int): Maximum number of retry attempts
        delay (int): Delay in seconds between retries
        
    Returns:
        requests.Response: The response object
        
    Raises:
        Exception: If all retries fail
    """
    for attempt in range(max_retries):
        try:
            response = requests.get(url)
            response.raise_for_status()
            return response
        except requests.exceptions.RequestException as e:
            if attempt < max_retries - 1:
                print(f"Request failed (attempt {attempt + 1}/{max_retries}): {e}")
                print(f"Retrying in {delay} seconds...")
                time.sleep(delay)
            else:
                raise e

def get_altitude(x, y):
    """
    Queries the IGN WMS elevation service for a given x, y coordinate in EPSG:2154.
    Uses tile caching to minimize network requests.

    Args:
        x (float): The x-coordinate in EPSG:2154.
        y (float): The y-coordinate in EPSG:2154.

    Returns:
        float or None: The altitude in meters, or None if the query fails.
    """
    try:
        # Calculate tile boundaries (10km x 10km)
        x_min = (x // TILE_SIZE) * TILE_SIZE
        y_min = (y // TILE_SIZE) * TILE_SIZE
        x_max = x_min + TILE_SIZE
        y_max = y_min + TILE_SIZE
        
        key = (x_min, y_min)
        
        if key not in tile_cache:
            # Check if tile is cached on disk
            if persistent_temp_dir:
                cache_file = os.path.join(persistent_temp_dir, 'tile_cache', f"{x_min}_{y_min}.npy")
                if os.path.exists(cache_file):
                    print(f"Loading tile from cache: {cache_file}")
                    data = np.load(cache_file)
                    tile_cache[key] = data
                else:
                    # Download and save
                    os.makedirs(os.path.dirname(cache_file), exist_ok=True)
                    bbox = f"{x_min},{y_min},{x_max},{y_max}"
                    print(f"Downloading tile for BBOX: {bbox}")
                    url = (
                        f"{WMS_URL}?SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap"
                        f"&BBOX={bbox}&CRS=EPSG:2154&WIDTH={TILE_PIXELS}&HEIGHT={TILE_PIXELS}"
                        f"&LAYERS={LAYER}&STYLES=&FORMAT=image/x-bil;bits=32"
                    )
                    
                    # Download the tile with retry
                    response = download_tile_with_retry(url)
                    
                    # Read as float32 array
                    data = np.frombuffer(response.content, dtype=np.float32).reshape((TILE_PIXELS, TILE_PIXELS))
                    np.save(cache_file, data)
                    tile_cache[key] = data
            else:
                # Download without saving
                bbox = f"{x_min},{y_min},{x_max},{y_max}"
                print(f"Downloading tile for BBOX: {bbox}")
                url = (
                    f"{WMS_URL}?SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap"
                    f"&BBOX={bbox}&CRS=EPSG:2154&WIDTH={TILE_PIXELS}&HEIGHT={TILE_PIXELS}"
                    f"&LAYERS={LAYER}&STYLES=&FORMAT=image/x-bil;bits=32"
                )
                
                # Download the tile with retry
                response = download_tile_with_retry(url)
                
                # Read as float32 array
                data = np.frombuffer(response.content, dtype=np.float32).reshape((TILE_PIXELS, TILE_PIXELS))
                tile_cache[key] = data
        
        # Get the data from cache
        data = tile_cache[key]
        
        # Calculate pixel coordinates within the tile
        pixel_x = int((x - x_min) / PIXEL_SIZE)
        pixel_y = int((y - y_min) / PIXEL_SIZE)
        
        # Clamp to valid range
        pixel_x = max(0, min(pixel_x, TILE_PIXELS - 1))
        pixel_y = max(0, min(pixel_y, TILE_PIXELS - 1))
        
        # Extract altitude value
        value = data[pixel_y, pixel_x]
        
        return float(value)
    
    except Exception as e:
        print(f"Error getting altitude for ({x}, {y}): {e}")
        return None

def process_geojson_file(file_path, dest_dir):
    """
    Processes a single GeoJSON file: adds altitude to each point feature and saves to destination.

    Args:
        file_path (str): Path to the source GeoJSON file.
        dest_dir (str): Path to the destination directory.
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)

    # Ensure it's a FeatureCollection
    if geojson_data.get('type') != 'FeatureCollection':
        print(f"Skipping {file_path}: Not a FeatureCollection")
        return

    features = geojson_data['features']
    for feature in features:
        geometry = feature.get('geometry', {})
        if geometry.get('type') == 'Point':
            coords = geometry['coordinates']
            # GeoJSON coordinates are [x, y] in EPSG:2154 (Lambert 93)
            x, y = coords[0], coords[1]
            altitude = get_altitude(x, y)
            if altitude is not None:
                feature['properties']['ALTITUDE'] = altitude
            else:
                feature['properties']['ALTITUDE'] = None  # Or skip adding if preferred
        else:
            print(f"Warning: Feature with non-point geometry in {file_path}")

    # Save the modified GeoJSON to the destination directory
    base_name = os.path.basename(file_path)
    dest_path = os.path.join(dest_dir, base_name)
    with open(dest_path, 'w', encoding='utf-8') as f:
        json.dump(geojson_data, f, indent=2, ensure_ascii=False)
    print(f"Processed and saved: {dest_path}")

def process_shapefile(file_path, dest_dir):
    """
    Processes a single ESRI Shapefile: adds altitude to each point feature and saves as Shapefile.

    Args:
        file_path (str): Path to the source Shapefile (.shp).
        dest_dir (str): Path to the destination directory.
    """
    try:
        # Read the shapefile
        gdf = gpd.read_file(file_path)
        
        # Filter only Point geometries
        points_gdf = gdf[gdf.geometry.type == 'Point'].copy()
        
        if len(points_gdf) == 0:
            print(f"No point geometries found in {file_path}, skipping.")
            return
        
        # Ensure CRS is EPSG:2154 (Lambert 93)
        if points_gdf.crs is None:
            print(f"Warning: No CRS found in {file_path}, assuming EPSG:2154")
            points_gdf.set_crs("EPSG:2154", inplace=True)
        elif points_gdf.crs != "EPSG:2154":
            points_gdf = points_gdf.to_crs("EPSG:2154")
        
        # Add altitude to each point
        altitudes = []
        for geom in points_gdf.geometry:
            x, y = geom.x, geom.y
            alt = get_altitude(x, y)
            altitudes.append(alt)
        
        points_gdf['ALTITUDE'] = altitudes
        
        # Save as Shapefile
        base_name = Path(file_path).stem
        dest_path = os.path.join(dest_dir, base_name)
        points_gdf.to_file(dest_path + '.shp')
        print(f"Processed and saved: {dest_path}.shp")
        
    except Exception as e:
        print(f"Error processing shapefile {file_path}: {e}")

def process_zip_file(file_path, dest_dir):
    """
    Processes a BRGM ZIP file: extracts shapefiles, processes each .shp file, and re-creates the ZIP.

    Args:
        file_path (str): Path to the ZIP file.
        dest_dir (str): Path to the destination directory.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Extract ZIP
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
            
            # Find all .shp files
            shp_files = list(Path(temp_dir).rglob('*.shp'))
            
            processed_any = False
            for shp_file in shp_files:
                print(f"Processing shapefile from ZIP: {shp_file.name}")
                # Process in place
                temp_dest_dir = str(shp_file.parent)
                process_shapefile_inplace(str(shp_file), temp_dest_dir)
                processed_any = True
            
            if processed_any:
                # Create new ZIP with modified files
                # Extract code from file name (e.g., tmpXXX_074.zip -> 074)
                file_stem = Path(file_path).stem
                code = file_stem.split('_')[-1]  # Get the last part after _
                dest_zip_path = os.path.join(dest_dir, f"GEO050K_HARM_{code}.zip")
                with zipfile.ZipFile(dest_zip_path, 'w', zipfile.ZIP_DEFLATED) as zip_out:
                    for root, dirs, files in os.walk(temp_dir):
                        for file in files:
                            file_path_in_zip = os.path.join(root, file)
                            arcname = os.path.relpath(file_path_in_zip, temp_dir)
                            zip_out.write(file_path_in_zip, arcname)
                print(f"Processed and saved: {dest_zip_path}")
            else:
                print(f"No shapefiles with points found in {file_path}, skipping ZIP creation.")
                
        except Exception as e:
            print(f"Error processing ZIP file {file_path}: {e}")

def process_shapefile_inplace(file_path, dest_dir):
    """
    Processes a shapefile and saves it in the same directory (for ZIP processing).
    """
    try:
        # Read the shapefile
        gdf = gpd.read_file(file_path)
        
        # Filter only Point geometries
        points_gdf = gdf[gdf.geometry.type == 'Point'].copy()
        
        if len(points_gdf) == 0:
            return  # No points, leave as is
        
        # Ensure CRS is EPSG:2154 (Lambert 93)
        if points_gdf.crs is None:
            points_gdf.set_crs("EPSG:2154", inplace=True)
        elif points_gdf.crs != "EPSG:2154":
            points_gdf = points_gdf.to_crs("EPSG:2154")
        
        # Add altitude to each point
        altitudes = []
        for geom in points_gdf.geometry:
            x, y = geom.x, geom.y
            alt = get_altitude(x, y)
            altitudes.append(alt)
        
        points_gdf['ALTITUDE'] = altitudes
        
        # Save in the same location
        base_name = Path(file_path).stem
        dest_path = os.path.join(dest_dir, base_name)
        points_gdf.to_file(dest_path + '.shp')
        
    except Exception as e:
        print(f"Error processing shapefile {file_path}: {e}")


def download_and_process_brgm(dest_dir, departments=None, keep_temp=False):
    """
    Download BRGM department ZIP files and process each downloaded ZIP into `dest_dir` preserving ZIP output.
    """
    base_url = "https://infoterre.brgm.fr/telechargements/BDCharm50/GEO050K_HARM_"

    if departments:
        codes = departments
    else:
        # Build department codes: 001..095 excluding 020, and include 02A/02B for Corsica
        codes = [f"{i:03d}" for i in range(1, 96) if i != 20]
        codes.insert(1, '02A')
        codes.insert(2, '02B')

    # Convert to set for easier processing
    codes_set = set(codes) if codes else set()
    processed_codes = set()

    # Handle special combined ZIP files
    # Île-de-France (75, 77, 78, 91, 92, 93, 94, 95)
    ile_de_france_deps = {'075', '077', '078', '091', '092', '093', '094', '095'}
    if ile_de_france_deps.issubset(codes_set):
        zip_filename = "GEO050K_HARM_075_077_078_091_092_093_094_095.zip"
        if keep_temp:
            tmp_path = os.path.join(tempfile.gettempdir(), zip_filename)
            if os.path.exists(tmp_path):
                print(f"Reusing existing Île-de-France ZIP from {tmp_path}")
                try:
                    process_zip_file(tmp_path, dest_dir)
                    processed_codes.update(ile_de_france_deps)
                except Exception as e:
                    print(f"Error processing Île-de-France ZIP: {e}")
            else:
                # Download and save
                url = f"{base_url}075_077_078_091_092_093_094_095.zip"
                print(f"Downloading BRGM ZIP for Île-de-France from {url}")
                try:
                    resp = requests.get(url, stream=True, timeout=30)
                    if resp.status_code != 200:
                        print(f"  -> Not available (status {resp.status_code}), skipping Île-de-France")
                    else:
                        with open(tmp_path, 'wb') as f:
                            for chunk in resp.iter_content(chunk_size=8192):
                                if chunk:
                                    f.write(chunk)
                        # Process the ZIP
                        process_zip_file(tmp_path, dest_dir)
                        processed_codes.update(ile_de_france_deps)
                except Exception as e:
                    print(f"Error downloading or processing Île-de-France: {e}")
        else:
            # Download without saving
            url = f"{base_url}075_077_078_091_092_093_094_095.zip"
            print(f"Downloading BRGM ZIP for Île-de-France from {url}")
            try:
                resp = requests.get(url, stream=True, timeout=30)
                if resp.status_code != 200:
                    print(f"  -> Not available (status {resp.status_code}), skipping Île-de-France")
                else:
                    with tempfile.NamedTemporaryFile(suffix="_ile_de_france.zip", delete=False) as tmp:
                        for chunk in resp.iter_content(chunk_size=8192):
                            if chunk:
                                tmp.write(chunk)
                        tmp_path = tmp.name
                    # Process the ZIP
                    process_zip_file(tmp_path, dest_dir)
                    processed_codes.update(ile_de_france_deps)
            except Exception as e:
                print(f"Error downloading or processing Île-de-France: {e}")
            finally:
                if 'tmp_path' in locals() and os.path.exists(tmp_path):
                    try:
                        os.unlink(tmp_path)
                    except Exception:
                        pass

    # Nord-Pas-de-Calais (59, 62)
    if '059' in codes_set and '062' in codes_set:
        zip_filename = "GEO050K_HARM_059_062.zip"
        if keep_temp:
            tmp_path = os.path.join(tempfile.gettempdir(), zip_filename)
            if os.path.exists(tmp_path):
                print(f"Reusing existing Nord-Pas-de-Calais ZIP from {tmp_path}")
                try:
                    process_zip_file(tmp_path, dest_dir)
                    processed_codes.add('059')
                    processed_codes.add('062')
                except Exception as e:
                    print(f"Error processing Nord-Pas-de-Calais ZIP: {e}")
            else:
                # Download and save
                url = f"{base_url}059_062.zip"
                print(f"Downloading BRGM ZIP for Nord-Pas-de-Calais from {url}")
                try:
                    resp = requests.get(url, stream=True, timeout=30)
                    if resp.status_code != 200:
                        print(f"  -> Not available (status {resp.status_code}), skipping Nord-Pas-de-Calais")
                    else:
                        with open(tmp_path, 'wb') as f:
                            for chunk in resp.iter_content(chunk_size=8192):
                                if chunk:
                                    f.write(chunk)
                        # Process the ZIP
                        process_zip_file(tmp_path, dest_dir)
                        processed_codes.add('059')
                        processed_codes.add('062')
                except Exception as e:
                    print(f"Error downloading or processing Nord-Pas-de-Calais: {e}")
        else:
            # Download without saving
            url = f"{base_url}059_062.zip"
            print(f"Downloading BRGM ZIP for Nord-Pas-de-Calais from {url}")
            try:
                resp = requests.get(url, stream=True, timeout=30)
                if resp.status_code != 200:
                    print(f"  -> Not available (status {resp.status_code}), skipping Nord-Pas-de-Calais")
                else:
                    with tempfile.NamedTemporaryFile(suffix="_nord_pas_de_calais.zip", delete=False) as tmp:
                        for chunk in resp.iter_content(chunk_size=8192):
                            if chunk:
                                tmp.write(chunk)
                        tmp_path = tmp.name
                    # Process the ZIP
                    process_zip_file(tmp_path, dest_dir)
                    processed_codes.add('059')
                    processed_codes.add('062')
            except Exception as e:
                print(f"Error downloading or processing Nord-Pas-de-Calais: {e}")
            finally:
                if 'tmp_path' in locals() and os.path.exists(tmp_path):
                    try:
                        os.unlink(tmp_path)
                    except Exception:
                        pass

    # Process remaining individual departments
    for code in codes:
        if code in processed_codes:
            continue  # Already processed as part of combined ZIP

        zip_filename = f"GEO050K_HARM_{code}.zip"
        if keep_temp:
            tmp_path = os.path.join(tempfile.gettempdir(), zip_filename)
            if os.path.exists(tmp_path):
                print(f"Reusing existing ZIP for department {code} from {tmp_path}")
                try:
                    process_zip_file(tmp_path, dest_dir)
                except Exception as e:
                    print(f"Error processing ZIP for department {code}: {e}")
                continue
        
        url = f"{base_url}{code}.zip"
        print(f"Downloading BRGM ZIP for department {code} from {url}")
        try:
            resp = requests.get(url, stream=True, timeout=30)
            if resp.status_code != 200:
                print(f"  -> Not available (status {resp.status_code}), skipping {code}")
                continue

            if keep_temp:
                tmp_path = os.path.join(tempfile.gettempdir(), zip_filename)
                with open(tmp_path, 'wb') as f:
                    for chunk in resp.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
            else:
                with tempfile.NamedTemporaryFile(suffix=f"_{code}.zip", delete=False) as tmp:
                    for chunk in resp.iter_content(chunk_size=8192):
                        if chunk:
                            tmp.write(chunk)
                    tmp_path = tmp.name

            # Process the ZIP and write resulting ZIP to dest_dir
            process_zip_file(tmp_path, dest_dir)

        except Exception as e:
            print(f"Error downloading or processing department {code}: {e}")
        finally:
            if not keep_temp:
                try:
                    if 'tmp_path' in locals() and os.path.exists(tmp_path):
                        os.unlink(tmp_path)
                except Exception:
                    pass
        zip_filename = f"GEO050K_HARM_{code}.zip"
        if keep_temp:
            tmp_path = os.path.join(tempfile.gettempdir(), zip_filename)
            if os.path.exists(tmp_path):
                print(f"Reusing existing ZIP for department {code} from {tmp_path}")
                try:
                    process_zip_file(tmp_path, dest_dir)
                except Exception as e:
                    print(f"Error processing ZIP for department {code}: {e}")
                continue
        
        url = f"{base_url}{code}.zip"
        print(f"Downloading BRGM ZIP for department {code} from {url}")
        try:
            resp = requests.get(url, stream=True, timeout=30)
            if resp.status_code != 200:
                print(f"  -> Not available (status {resp.status_code}), skipping {code}")
                continue

            if keep_temp:
                tmp_path = os.path.join(tempfile.gettempdir(), zip_filename)
                with open(tmp_path, 'wb') as f:
                    for chunk in resp.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
            else:
                with tempfile.NamedTemporaryFile(suffix=f"_{code}.zip", delete=False) as tmp:
                    for chunk in resp.iter_content(chunk_size=8192):
                        if chunk:
                            tmp.write(chunk)
                    tmp_path = tmp.name

            # Process the ZIP and write resulting ZIP to dest_dir
            process_zip_file(tmp_path, dest_dir)

        except Exception as e:
            print(f"Error downloading or processing department {code}: {e}")
        finally:
            if not keep_temp:
                try:
                    if 'tmp_path' in locals() and os.path.exists(tmp_path):
                        os.unlink(tmp_path)
                except Exception:
                    pass



def main():
    """Main CLI: supports processing a source directory or downloading all BRGM departments with --from-brgm.

    Usage patterns:
      python add_altitude.py --in <source_dir> --out <dest_dir>
      python add_altitude.py --from-brgm --out <dest_dir> [--departments <codes>]
    """
    import argparse

    parser = argparse.ArgumentParser(description='Add altitude to BRGM data (GeoJSON, Shapefile, ZIP).')
    parser.add_argument('--in', dest='source', help='Source directory to process')
    parser.add_argument('--out', dest='dest', required=True, help='Destination directory')
    parser.add_argument('--from-brgm', action='store_true', help='Download BRGM department ZIPs and process them')
    parser.add_argument('--departments', help='Department codes to download (comma-separated, or "all" for all departments)')
    parser.add_argument('--temp-dir', help='Temporary directory to use (will be preserved if specified)')

    args = parser.parse_args()

    # Set temp dir
    if args.temp_dir:
        global persistent_temp_dir
        persistent_temp_dir = args.temp_dir
        os.makedirs(persistent_temp_dir, exist_ok=True)
        # Use persistent temp dir
        import tempfile
        tempfile.tempdir = persistent_temp_dir
    else:
        temp_dir = None

    if args.from_brgm:
        if not args.dest:
            print('Error: --out is required')
            sys.exit(1)
        dest_dir = args.dest
        os.makedirs(dest_dir, exist_ok=True)
        
        # Parse departments
        if args.departments:
            if args.departments.lower() == 'all':
                departments = None  # All
            else:
                departments = [code.strip() for code in args.departments.split(',')]
        else:
            departments = None  # All
        
        download_and_process_brgm(dest_dir, departments, keep_temp=bool(args.temp_dir))
    else:
        # Process local files
        if not args.source or not args.dest:
            print('Usage: python add_altitude.py --in <source_dir> --out <dest_dir>')
            sys.exit(1)
        
        source_dir = args.source
        dest_dir = args.dest
        
        # Check if source directory exists
        if not os.path.isdir(source_dir):
            print(f"Error: Source directory '{source_dir}' does not exist.")
            sys.exit(1)
        
        # Create destination directory if it doesn't exist
        os.makedirs(dest_dir, exist_ok=True)
        
        # Process each supported file in the source directory
        for filename in os.listdir(source_dir):
            file_path = os.path.join(source_dir, filename)
            if filename.endswith('.geojson'):
                process_geojson_file(file_path, dest_dir)
            elif filename.endswith('.shp'):
                process_shapefile(file_path, dest_dir)
            elif filename.endswith('.zip'):
                process_zip_file(file_path, dest_dir)
            else:
                print(f"Skipping unsupported file: {filename}")

if __name__ == "__main__":
    main()