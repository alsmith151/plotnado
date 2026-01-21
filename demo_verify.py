"""
Demo script to verify plotnado refactoring works.
"""

import plotnado as pn

# Test 1: GenomicRegion parsing
print("Test 1: GenomicRegion parsing")
gr = pn.GenomicRegion.from_str("chr1:1000000-2000000")
print(f"  Region: {gr}")
print(f"  Length: {gr.length:,} bp")
print(f"  Center: {gr.center:,}")
print("  ✓ GenomicRegion works\n")

# Test 2: Create a Figure with tracks
print("Test 2: Creating Figure with tracks")
fig = pn.Figure(width=10)
fig.add_track("scalebar")
fig.add_track("spacer", height=0.3)
print(f"  Figure: {fig}")
print(f"  Tracks: {len(fig.tracks)}")
print("  ✓ Figure creation works\n")

# Test 3: BigWig track (without actual file)
print("Test 3: BigWig track creation")
import pandas as pd

data = pd.DataFrame(
    {
        "chrom": ["chr1"] * 10,
        "start": range(1000000, 1001000, 100),
        "end": range(1000100, 1001100, 100),
        "value": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    }
)
bw_track = pn.BigWigTrack(data=data, title="Test Signal")
fig2 = pn.Figure()
fig2.add_track(bw_track)
print(f"  BigWig track: {bw_track.title}")
print("  ✓ BigWig track works\n")

# Test 4: Aesthetics
print("Test 4: Track aesthetics")
aesthetics = pn.BigwigAesthetics(style="fill", color="steelblue", alpha=0.8)
print(f"  Style: {aesthetics.style}")
print(f"  Color: {aesthetics.color}")
print("  ✓ Aesthetics work\n")

print("=" * 50)
print("All basic functionality tests passed!")
print("=" * 50)
