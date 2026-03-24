#!/usr/bin/env python3
"""
Grape_DRF Loader (ZIP OR pre-extracted folder support)

Author: N6RFM + ChatGPT

PURPOSE
-------
This script prepares GRAPE DRF datasets for use in the G3ZIL grapeDRF_doppler_model
processing pipeline.

It supports TWO input types:
    1) Raw OBS*.zip files
    2) Already-extracted dataset folders (ch0_*)

The script ensures:
    - ZIPs are flattened correctly (no extra "ch0/" nesting)
    - Existing folders are reused without modification
    - Config file timestamp is updated
    - Working directory is refreshed cleanly

PIPELINE OVERVIEW
-----------------
1. Scan library directory for:
       - OBS*.zip
       - ch0_* folders
2. User selects input source
3. User enters callsign
4. If ZIP:
       - Extract into new ch0_<CALLSIGN>_<TIMESTAMP> folder
       - Strip top-level directory from ZIP
5. If folder:
       - Use folder directly
6. Update config/<CALLSIGN>_config.ini
7. Clear working directory
8. Copy dataset into working directory
"""

import os
import shutil
import re
import zipfile

# -----------------------------
# PATH CONFIGURATION
# -----------------------------
# Base directory of your Doppler model project
BASE_DIR = "/home/bob/grapeDRF_doppler_model"

# Where ZIP files and extracted datasets live
LIBRARY_DIR = os.path.join(BASE_DIR, "data/library")

# Working directory used by processing pipeline
DATA_ROOT = os.path.join(BASE_DIR, "data/psws_grapeDRF")

# Config files directory
CONFIG_DIR = os.path.join(BASE_DIR, "config")


# -----------------------------
# FIND INPUT SOURCES
# -----------------------------
def find_sources():
    """
    Scan LIBRARY_DIR and return all valid input sources.

    Returns:
        List of tuples:
            ("zip", filename)
            ("folder", foldername)

    Logic:
        - ZIP files must match OBS*.zip
        - Folders must start with ch0_
    """
    sources = []

    for item in os.listdir(LIBRARY_DIR):
        full_path = os.path.join(LIBRARY_DIR, item)

        # Identify ZIP files
        if item.startswith("OBS") and item.endswith(".zip"):
            sources.append(("zip", item))

        # Identify already-extracted dataset folders
        elif os.path.isdir(full_path) and item.startswith("ch0_"):
            sources.append(("folder", item))

    return sources


# -----------------------------
# USER SELECTION
# -----------------------------
def select_source(sources):
    """
    Present available sources and allow user selection.

    Input:
        sources = list from find_sources()

    Returns:
        (type, name)
    """
    print("\nAvailable sources:\n")

    for i, (stype, name) in enumerate(sources):
        print(f"{i}: [{stype}] {name}")

    idx = int(input("\nSelect index: "))
    return sources[idx]


# -----------------------------
# ZIP TIMESTAMP PARSING
# -----------------------------
def parse_zip_datetime(zip_name):
    """
    Extract timestamp from ZIP filename.

    Expected format:
        OBSYYYY-MM-DDTHH-MM*.zip

    Example:
        OBS2026-03-20T12-30.zip

    Returns:
        (year, month, day, hour, minute) as integers
    """
    match = re.search(r"OBS(\d{4})-(\d{2})-(\d{2})T(\d{2})-(\d{2})", zip_name)

    if not match:
        raise ValueError("Invalid ZIP filename format")

    return tuple(map(int, match.groups()))


def format_dt(year, month, day, hour, minute):
    """
    Convert datetime components to compact string.

    Example:
        (2026,3,20,12,30) → "202603201230"
    """
    return f"{year:04d}{month:02d}{day:02d}{hour:02d}{minute:02d}"


# -----------------------------
# CALLSIGN INPUT
# -----------------------------
def get_callsign():
    """
    Prompt user for callsign.

    Default:
        N6RFM

    Returns:
        Uppercase callsign string
    """
    callsign = input("\nEnter callsign (default N6RFM): ").strip().upper()
    return callsign if callsign else "N6RFM"


# -----------------------------
# ZIP EXTRACTION (FLATTENED)
# -----------------------------
def extract_zip_flat(zip_file, callsign):
    """
    Extract ZIP into a clean dataset folder.

    CRITICAL BEHAVIOR:
        - Removes top-level directory (e.g., "ch0/")
        - Preserves all internal structure below that
        - Prevents unwanted nested "ch0/ch0/" issues

    Output folder:
        ch0_<CALLSIGN>_<YYYYMMDDHHMM>

    Returns:
        (folder_name, (year, month, day, hour, minute))
    """

    # Extract timestamp from filename
    year, month, day, hour, minute = parse_zip_datetime(zip_file)
    dt_str = format_dt(year, month, day, hour, minute)

    folder_name = f"ch0_{callsign}_{dt_str}"
    dest_path = os.path.join(LIBRARY_DIR, folder_name)

    os.makedirs(dest_path, exist_ok=True)

    zip_path = os.path.join(LIBRARY_DIR, zip_file)

    print(f"\nExtracting ZIP → {folder_name}")

    with zipfile.ZipFile(zip_path, 'r') as zip_ref:

        # ----------------------------------------
        # Detect if ZIP has a single top-level folder
        # Example:
        #   ch0/file1
        #   ch0/file2
        #
        # If so, we remove it
        # ----------------------------------------
        top_levels = set()

        for m in zip_ref.infolist():
            if m.filename.strip():
                top_levels.add(m.filename.split('/')[0])

        top_level = list(top_levels)[0] if len(top_levels) == 1 else None

        # ----------------------------------------
        # Extract files manually
        # ----------------------------------------
        for m in zip_ref.infolist():

            # Skip directory entries
            if m.is_dir():
                continue

            parts = m.filename.split('/')

            # Remove top-level folder if present
            if top_level and parts[0] == top_level:
                parts = parts[1:]

            # Skip empty paths
            if not parts:
                continue

            target = os.path.join(dest_path, *parts)

            # Ensure parent directories exist
            os.makedirs(os.path.dirname(target), exist_ok=True)

            # Copy file from ZIP → disk
            with zip_ref.open(m) as src, open(target, 'wb') as dst:
                shutil.copyfileobj(src, dst)

    print("Extraction complete (flattened)")
    return folder_name, (year, month, day, hour, minute)


# -----------------------------
# CONFIG UPDATE
# -----------------------------
def update_config(callsign, year, month, day, hour, minute):
    """
    Update UT timestamp inside config file.

    File:
        config/<CALLSIGN>_config.ini

    Replaces line starting with:
        ut = [...]

    With:
        ut = [year,month,day,hour,minute]
    """
    config_file = os.path.join(CONFIG_DIR, f"{callsign}_config.ini")

    if not os.path.exists(config_file):
        raise FileNotFoundError(config_file)

    new_ut = f"ut = [{year},{month},{day},{hour},{minute}]"

    with open(config_file) as f:
        lines = f.readlines()

    with open(config_file, "w") as f:
        for line in lines:
            if line.strip().startswith("ut"):
                f.write(new_ut + "\n")
            else:
                f.write(line)

    print(f"Updated config: {config_file}")


# -----------------------------
# CLEAR WORKING DIRECTORY
# -----------------------------
def clear_destination(dest_dir):
    """
    Delete all contents of working directory AFTER user confirmation.

    This prevents accidental data loss.
    """
    confirm = input(f"\n⚠️ Clear {dest_dir}? (y/n): ").lower()

    if confirm != "y":
        print("Aborted")
        exit(0)

    for item in os.listdir(dest_dir):
        path = os.path.join(dest_dir, item)

        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.unlink(path)


# -----------------------------
# COPY DATA TO WORKING DIRECTORY
# -----------------------------
def copy_data(folder_name, callsign):
    """
    Copy dataset into working directory.

    Source:
        data/library/<folder_name>

    Destination:
        data/psws_grapeDRF/ch0_<CALLSIGN>

    Behavior:
        - Clears destination first
        - Copies everything (no extra nesting)
    """
    src = os.path.join(LIBRARY_DIR, folder_name)
    dest = os.path.join(DATA_ROOT, f"ch0_{callsign}")

    os.makedirs(dest, exist_ok=True)

    clear_destination(dest)

    print(f"\nCopying → {dest}")

    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dest, item)

        if os.path.isdir(s):
            shutil.copytree(s, d)
        else:
            shutil.copy2(s, d)

    print("Copy complete")
    return dest


# -----------------------------
# MAIN PROGRAM
# -----------------------------
def main():
    """
    Main control flow.

    Handles:
        - Source selection
        - ZIP vs folder branching
        - Config update
        - Data deployment
    """
    try:
        sources = find_sources()

        if not sources:
            print("No ZIPs or folders found")
            return

        stype, name = select_source(sources)
        callsign = get_callsign()

        # ----------------------------------------
        # ZIP INPUT PATH
        # ----------------------------------------
        if stype == "zip":
            folder_name, (y, m, d, h, mi) = extract_zip_flat(name, callsign)

        # ----------------------------------------
        # FOLDER INPUT PATH
        # ----------------------------------------
        else:
            folder_name = name

            # Extract timestamp from folder name
            # Expected format:
            #   ch0_CALLSIGN_YYYYMMDDHHMM
            match = re.search(r"(\d{12})$", name)

            if not match:
                raise ValueError("Folder missing timestamp")

            dt = match.group(1)

            y = int(dt[0:4])
            m = int(dt[4:6])
            d = int(dt[6:8])
            h = int(dt[8:10])
            mi = int(dt[10:12])

        # ----------------------------------------
        # UPDATE CONFIG
        # ----------------------------------------
        update_config(callsign, y, m, d, h, mi)

        # ----------------------------------------
        # COPY DATA INTO WORKING AREA
        # ----------------------------------------
        dest = copy_data(folder_name, callsign)

        print(f"\n✅ READY → {dest}")

    except Exception as e:
        print(f"ERROR: {e}")


# -----------------------------
# ENTRY POINT
# -----------------------------
if __name__ == "__main__":
    main()
