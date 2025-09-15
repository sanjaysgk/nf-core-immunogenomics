#!/usr/bin/env python3
"""
Download test data for nf-core/immunogenomics pipeline
Texas Cancer Research Biobank (TCRB) dataset from nf-core/test-datasets

Usage:
    python bin/download_test_data.py --outdir assets/test_data [--data-type wes|rnaseq|all]
"""

import argparse
import os
import sys
import urllib.request
from pathlib import Path

# Test data URLs from nf-core/test-datasets
BASE_URL = "https://github.com/nf-core/test-datasets/raw/rnadnavar/data/tcrb"

TEST_FILES = {
    "wes": {
        "TCRBOA7-N-WEX.read1.250000.fastq.gz": f"{BASE_URL}/fastq/TCRBOA7-N-WEX.read1.250000.fastq.gz",
        "TCRBOA7-N-WEX.read2.250000.fastq.gz": f"{BASE_URL}/fastq/TCRBOA7-N-WEX.read2.250000.fastq.gz",
        "TCRBOA7-T-WEX.read1.250000.fastq.gz": f"{BASE_URL}/fastq/TCRBOA7-T-WEX.read1.250000.fastq.gz",
        "TCRBOA7-T-WEX.read2.250000.fastq.gz": f"{BASE_URL}/fastq/TCRBOA7-T-WEX.read2.250000.fastq.gz",
    },
    "rnaseq": {
        "TCRBOA7-T-RNA.read1.250000.fastq.gz": f"{BASE_URL}/fastq/TCRBOA7-T-RNA.read1.250000.fastq.gz",
        "TCRBOA7-T-RNA.read2.250000.fastq.gz": f"{BASE_URL}/fastq/TCRBOA7-T-RNA.read2.250000.fastq.gz",
    }
}

def download_file(url, output_path):
    """Download a file from URL to output path with progress indication"""
    try:
        print(f"Downloading {os.path.basename(output_path)}...")
        
        def progress_hook(block_num, block_size, total_size):
            if total_size > 0:
                percent = min(100, (block_num * block_size * 100) // total_size)
                print(f"\r  Progress: {percent}%", end='', flush=True)
        
        urllib.request.urlretrieve(url, output_path, progress_hook)
        print()  # New line after progress
        return True
    except Exception as e:
        print(f"\nError downloading {url}: {e}")
        return False

def create_samplesheets(output_dir):
    """Create samplesheets for the downloaded test data"""
    
    # WES samplesheet
    wes_samplesheet = output_dir / "samplesheet_wes.csv"
    with open(wes_samplesheet, 'w') as f:
        f.write("sample,fastq_1,fastq_2,sample_type,patient_id\n")
        f.write(f"TCRBOA7_normal,{output_dir}/fastq/TCRBOA7-N-WEX.read1.250000.fastq.gz,{output_dir}/fastq/TCRBOA7-N-WEX.read2.250000.fastq.gz,normal,TCRBOA7\n")
        f.write(f"TCRBOA7_tumor,{output_dir}/fastq/TCRBOA7-T-WEX.read1.250000.fastq.gz,{output_dir}/fastq/TCRBOA7-T-WEX.read2.250000.fastq.gz,tumor,TCRBOA7\n")
    
    # RNA-seq samplesheet
    rnaseq_samplesheet = output_dir / "samplesheet_rnaseq.csv"
    with open(rnaseq_samplesheet, 'w') as f:
        f.write("sample,fastq_1,fastq_2,condition,batch\n")
        f.write(f"TCRBOA7_tumor_rna,{output_dir}/fastq/TCRBOA7-T-RNA.read1.250000.fastq.gz,{output_dir}/fastq/TCRBOA7-T-RNA.read2.250000.fastq.gz,tumor,1\n")
    
    return wes_samplesheet, rnaseq_samplesheet

def main():
    parser = argparse.ArgumentParser(
        description="Download test data for nf-core/immunogenomics pipeline"
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default="assets/test_data",
        help="Output directory for test data (default: assets/test_data)"
    )
    parser.add_argument(
        "--data-type",
        choices=["wes", "rnaseq", "all"],
        default="all",
        help="Type of data to download (default: all)"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing files"
    )
    
    args = parser.parse_args()
    
    # Create output directories
    fastq_dir = args.outdir / "fastq"
    fastq_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 80)
    print("nf-core/immunogenomics Test Data Downloader")
    print("=" * 80)
    print("\nDataset: Texas Cancer Research Biobank (TCRB)")
    print("Source: nf-core/test-datasets")
    print("⚠️  This data is for testing purposes only")
    print("⚠️  Do not attempt to re-identify participants")
    print()
    
    # Determine which files to download
    files_to_download = {}
    if args.data_type in ["wes", "all"]:
        files_to_download.update(TEST_FILES["wes"])
    if args.data_type in ["rnaseq", "all"]:
        files_to_download.update(TEST_FILES["rnaseq"])
    
    # Download files
    success_count = 0
    total_files = len(files_to_download)
    
    for filename, url in files_to_download.items():
        output_path = fastq_dir / filename
        
        # Skip if file exists and not forcing
        if output_path.exists() and not args.force:
            print(f"Skipping {filename} (already exists, use --force to overwrite)")
            success_count += 1
            continue
        
        if download_file(url, output_path):
            success_count += 1
        else:
            print(f"Failed to download {filename}")
    
    # Create samplesheets
    if success_count > 0:
        print(f"\nCreating samplesheets...")
        wes_sheet, rnaseq_sheet = create_samplesheets(args.outdir)
        print(f"Created: {wes_sheet}")
        print(f"Created: {rnaseq_sheet}")
    
    # Summary
    print(f"\n{'='*50}")
    print(f"Download Summary:")
    print(f"  Downloaded: {success_count}/{total_files} files")
    print(f"  Output directory: {args.outdir}")
    
    if success_count == total_files:
        print("\n✅ Test data download completed successfully!")
        print("\nTo test the pipeline:")
        print(f"  WES:     nextflow run . --input {args.outdir}/samplesheet_wes.csv")
        print(f"  RNA-seq: nextflow run . --input {args.outdir}/samplesheet_rnaseq.csv")
    else:
        print(f"\n❌ Some downloads failed. Check network connection and try again.")
        sys.exit(1)

if __name__ == "__main__":
    main()