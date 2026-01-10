#!/usr/bin/env python3
"""
Diagnostic script for troubleshooting FeaVar alignment loading issues.

Usage:
    python diagnose.py data/sample_alignment.clw
"""

import sys
from pathlib import Path


def diagnose_alignment(filepath: str):
    """Run diagnostics on an alignment file."""
    
    print(f"=" * 60)
    print("FeaVar Alignment Diagnostic Tool")
    print(f"=" * 60)
    
    path = Path(filepath)
    
    # Check file exists
    print(f"\n1. Checking file exists...")
    if not path.exists():
        print(f"   ERROR: File not found: {filepath}")
        return
    print(f"   OK: File exists at {path.absolute()}")
    
    # Check file size
    print(f"\n2. Checking file size...")
    size = path.stat().st_size
    print(f"   File size: {size} bytes")
    if size == 0:
        print("   ERROR: File is empty!")
        return
    
    # Check file content
    print(f"\n3. Checking file content...")
    try:
        content = path.read_text()
        lines = content.split('\n')
        print(f"   Number of lines: {len(lines)}")
        print(f"   First line: {lines[0][:60]}...")
    except Exception as e:
        print(f"   ERROR reading file: {e}")
        return
    
    # Detect format
    print(f"\n4. Detecting format...")
    first_line = lines[0].strip()
    if first_line.startswith('>'):
        detected_format = 'fasta'
    elif 'CLUSTAL' in first_line.upper():
        detected_format = 'clustal'
    elif first_line.startswith('#'):
        detected_format = 'stockholm'
    else:
        detected_format = 'unknown'
    print(f"   Detected format: {detected_format}")
    
    # Try loading with BioPython
    print(f"\n5. Testing BioPython loading...")
    try:
        from Bio import AlignIO
        print(f"   BioPython imported successfully")
        
        # Try detected format
        formats_to_try = [detected_format] if detected_format != 'unknown' else []
        formats_to_try.extend(['clustal', 'fasta', 'phylip'])
        
        for fmt in formats_to_try:
            try:
                print(f"\n   Trying format: {fmt}")
                alignment = AlignIO.read(filepath, fmt)
                print(f"   SUCCESS with format '{fmt}'!")
                print(f"   Number of sequences: {len(alignment)}")
                print(f"   Alignment length: {alignment.get_alignment_length()}")
                print(f"   Sequence IDs: {[r.id for r in alignment]}")
                return fmt
            except Exception as e:
                error_str = str(e)
                if not error_str:
                    error_str = f"{type(e).__name__} (no message)"
                print(f"   Failed with {fmt}: {error_str}")
        
        print("\n   ERROR: Could not load with any format")
        
    except ImportError:
        print("   ERROR: BioPython not installed!")
        print("   Install with: pip install biopython")
    
    # Show raw file content for debugging
    print(f"\n6. Raw file content (first 20 lines):")
    print("-" * 40)
    for i, line in enumerate(lines[:20]):
        # Show special characters
        display = line.replace('\t', '→').replace(' ', '·')
        print(f"   {i+1:3}: {display[:70]}")
    print("-" * 40)
    
    return None


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python diagnose.py <alignment_file>")
        print("Example: python diagnose.py data/sample_alignment.clw")
        sys.exit(1)
    
    diagnose_alignment(sys.argv[1])
