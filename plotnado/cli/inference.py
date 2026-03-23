"""
Inference engine for track type detection and template generation heuristics.

This module provides:
- Track type classification from file paths and URLs
- Title inference from filenames
- Grouping heuristics for replicates and shared analysis
- Seqnado pipeline pattern detection (SAMPLE_ANTIBODY.bw)
- Confidence scoring for inference decisions
"""

import re
from pathlib import Path
from typing import Optional

from plotnado.template import TemplateTrackType


class SeqnadoPattern:
    """Detection and parsing for seqnado pipeline outputs."""
    
    # Standard seqnado output: SAMPLE-NAME_ANTIBODY.bigwig (without explicit bw/bigwig in name)
    # Examples: control_H3K27ac.bw, sample1_Input.bigwig
    # Won't match: THP1H3K4me1_bigWig.bigWig (has bigwig in filename before extension)
    PATTERN = re.compile(
        r'([A-Za-z0-9\-_]+)_([A-Za-z0-9]+)\.(?:bw|bigwig)$',
        re.IGNORECASE
    )
    
    # File type suffixes to exclude from antibody name
    EXCLUDED_ANTIBODIES = {'bigwig', 'bw', 'bigbed', 'bed', 'sorted', 'bam'}
    
    @staticmethod
    def is_seqnado(filename: str) -> bool:
        """Check if a filename matches seqnado pattern."""
        match = SeqnadoPattern.PATTERN.match(Path(filename).name)
        if not match:
            return False
        # Reject if antibody looks like a file type suffix
        _, antibody = match.groups()
        return antibody.lower() not in SeqnadoPattern.EXCLUDED_ANTIBODIES
    
    @staticmethod
    def parse(filename: str) -> Optional[tuple[str, str]]:
        """
        Parse a seqnado filename into (sample_name, antibody).
        
        Returns:
            Tuple of (sample_name, antibody) or None if not a seqnado file
        """
        match = SeqnadoPattern.PATTERN.search(Path(filename).name)
        if not match:
            return None
        
        sample, antibody = match.groups()
        
        # Reject if antibody is a file type
        if antibody.lower() in SeqnadoPattern.EXCLUDED_ANTIBODIES:
            return None
        
        return sample, antibody


class TrackClassifier:
    """Classifies data sources into track types."""
    
    # File extension to track type mappings
    EXTENSION_MAP = {
        '.bw': TemplateTrackType.BIGWIG,
        '.bigwig': TemplateTrackType.BIGWIG,
        '.bedgraph': TemplateTrackType.BEDGRAPH,
        '.bed': TemplateTrackType.BED,
        '.bigbed': TemplateTrackType.BED,  # BigBed is a BED format variant
        '.narrowpeak': TemplateTrackType.NARROWPEAK,
        '.broadpeak': TemplateTrackType.NARROWPEAK,  # Similar to narrowpeak
        '.bedpe': TemplateTrackType.LINKS,
        '.links': TemplateTrackType.LINKS,
    }
    
    # URL pattern mappings
    URL_PATTERNS = {
        r'\.bw$': TemplateTrackType.BIGWIG,
        r'\.bigwig$': TemplateTrackType.BIGWIG,
        r'\.bigbed$': TemplateTrackType.BED,  # BigBed format
        r'\.bed': TemplateTrackType.BED,
        r'\.narrowpeak': TemplateTrackType.NARROWPEAK,
        r'\.bedgraph': TemplateTrackType.BEDGRAPH,
    }
    
    @staticmethod
    def classify(path: str) -> tuple[TemplateTrackType, float]:
        """
        Classify a data source into a track type.
        
        Returns:
            Tuple of (track_type, confidence) where confidence is 0-1
        """
        # Check for URLs
        if path.startswith(('http://', 'https://', 's3://', 'ftp://')):
            for pattern, track_type in TrackClassifier.URL_PATTERNS.items():
                if re.search(pattern, path, re.IGNORECASE):
                    return track_type, 0.9  # High confidence for URL patterns
            return TemplateTrackType.UNKNOWN, 0.1
        
        # Check file extension
        path_lower = path.lower()
        for ext, track_type in TrackClassifier.EXTENSION_MAP.items():
            if path_lower.endswith(ext):
                return track_type, 0.95  # Very high confidence for direct extension match
        
        # Check Path object extension
        try:
            p = Path(path)
            suffix = p.suffix.lower()
            if suffix in TrackClassifier.EXTENSION_MAP:
                return TrackClassifier.EXTENSION_MAP[suffix], 0.95
        except (TypeError, ValueError):
            pass
        
        return TemplateTrackType.UNKNOWN, 0.0


class TitleInference:
    """Infers meaningful titles from file paths and filenames."""
    
    @staticmethod
    def infer(path: str, track_type: Optional[TemplateTrackType] = None) -> tuple[str, bool]:
        """
        Infer a title from a file path or URL.
        
        Args:
            path: File path or URL
            track_type: Optional track type to add appropriate suffix
        
        Returns:
            Tuple of (title, was_inferred)
        """
        try:
            # Extract filename from path or URL
            if '://' in path:  # URL
                path = path.split('/')[-1]
            
            filename = Path(path).stem  # Remove extension
            
            # Check for seqnado pattern first
            seqnado = SeqnadoPattern.parse(Path(path).name)
            if seqnado:
                sample, antibody = seqnado
                # Format as "Antibody (Sample)" for clarity
                title = f"{antibody} ({sample})"
                return title, True
            
            # Clean up common filename patterns
            title = TitleInference._clean_filename(filename)
            
            # Add track type suffix for clarity
            if title and track_type == TemplateTrackType.BED:
                title = f"{title} peaks"
            
            if title:
                return title, True
            return filename or "Unknown", False
        except (TypeError, ValueError):
            return "Unknown", False
    
    @staticmethod
    def _clean_filename(name: str) -> str:
        """
        Clean up common filename patterns while preserving scientific notation.
        
        Examples:
            "THP1H3K4me3_bigBed" → "THP1 H3K4me3"
            "control_rep1_Input" → "Control Input"
            "sample1_H3K27ac_sorted" → "Sample1 H3K27ac"
        """
        # Remove common file type suffixes (bigWig, bigBed, sorted, etc.)
        name = re.sub(r'[-_](bigwig|bigbed|sorted|bam|fastq|fq|txt)$', '', name, flags=re.IGNORECASE)
        
        # Separate camelCase words (e.g., "THP1H3K4me3" → "THP1 H3K4me3")
        # Insert space before uppercase letter that follows lowercase
        name = re.sub(r'([a-z])([A-Z])', r'\1 \2', name)
        
        # Insert space before histone marks (H3K, H4K, etc.)
        name = re.sub(r'([A-Za-z0-9])([H][0-9])', r'\1 \2', name)
        
        # Remove leading numbers and common control keywords
        name = re.sub(r'^[0-9]+[-_.]', '', name)
        name = re.sub(r'\b(rep|replicate|sample|input|control)[-_]?', '', name, flags=re.IGNORECASE)
        
        # Replace remaining underscores and hyphens with spaces
        name = re.sub(r'[_-]+', ' ', name)
        
        # Clean up multiple spaces
        name = re.sub(r'\s+', ' ', name)
        
        # Capitalize first letter of each word, but preserve existing capitals
        # (e.g., "H3K4me3" stays, "control" → "Control")
        words = name.split()
        capitalized = []
        for word in words:
            # If word starts with capital, keep original case
            if word and word[0].isupper():
                capitalized.append(word)
            else:
                # Otherwise capitalize first letter
                capitalized.append(word.capitalize() if word else '')
        
        name = ' '.join(capitalized)
        
        return name.strip()


class GroupingHeuristic:
    """Infers grouping for shared autoscaling and coloring."""
    
    @staticmethod
    def group_by_patterns(paths: list[str]) -> dict[str, list[int]]:
        """
        Group track indices by filename patterns.
        
        Looks for seqnado patterns first, then common prefixes, tokens, and replication patterns.
        Returns mapping of group_id -> list of track indices.
        """
        if not paths:
            return {}
        
        # Check if all are seqnado files
        seqnado_files = [SeqnadoPattern.parse(Path(p).name) for p in paths]
        all_seqnado = all(s is not None for s in seqnado_files)
        
        groups: dict[str, list[int]] = {}
        
        if all_seqnado:
            # Group by sample name (each antibody for the same sample shares scaling)
            samples: dict[str, list[int]] = {}
            for i, (sample, antibody) in enumerate(seqnado_files):
                if sample not in samples:
                    samples[sample] = []
                samples[sample].append(i)
            
            # Create groups for samples with multiple antibodies
            for sample, indices in samples.items():
                if len(indices) > 1:
                    group_id = f"{sample}_autoscale"
                    groups[group_id] = sorted(indices)
            
            return groups
        
        # Fall back to generic pattern matching
        stems = [Path(p).stem.lower() for p in paths]
        
        # Look for common patterns
        for i, stem1 in enumerate(stems):
            # Check if stem appears to be a replicate
            matches = []
            for j, stem2 in enumerate(stems):
                if i != j and GroupingHeuristic._are_replicates(stem1, stem2):
                    matches.append(j)
            
            if matches:
                # Found replicates - create a group
                group_id = GroupingHeuristic._extract_group_id(stem1)
                if group_id:
                    group_members = sorted([i] + matches)
                    group_key = f"{group_id}_group"
                    if group_key not in groups:
                        groups[group_key] = group_members
        
        return groups
    
    @staticmethod
    def group_by_sample(paths: list[str]) -> Optional[dict[str, list[int]]]:
        """
        Group by sample name (seqnado pattern only).
        
        Returns mapping of sample_name -> list of track indices,
        or None if files don't match seqnado pattern.
        """
        seqnado_files = [SeqnadoPattern.parse(Path(p).name) for p in paths]
        all_seqnado = all(s is not None for s in seqnado_files)
        
        if not all_seqnado:
            return None
        
        samples: dict[str, list[int]] = {}
        for i, (sample, antibody) in enumerate(seqnado_files):
            if sample not in samples:
                samples[sample] = []
            samples[sample].append(i)
        
        # Only return groups with multiple members (worth grouping)
        return {k: v for k, v in samples.items() if len(v) > 1}
    
    @staticmethod
    def group_by_antibody(paths: list[str]) -> Optional[dict[str, list[int]]]:
        """
        Group by antibody (seqnado pattern only).
        
        Returns mapping of antibody_name -> list of track indices,
        or None if files don't match seqnado pattern.
        """
        seqnado_files = [SeqnadoPattern.parse(Path(p).name) for p in paths]
        all_seqnado = all(s is not None for s in seqnado_files)
        
        if not all_seqnado:
            return None
        
        antibodies: dict[str, list[int]] = {}
        for i, (sample, antibody) in enumerate(seqnado_files):
            if antibody not in antibodies:
                antibodies[antibody] = []
            antibodies[antibody].append(i)
        
        # Only return groups with multiple members (worth grouping)
        return {k: v for k, v in antibodies.items() if len(v) > 1}
    
    @staticmethod
    def _are_replicates(stem1: str, stem2: str) -> bool:
        """Check if two stems appear to be replicates."""
        # Remove numeric suffixes (rep1, rep2, r1, r2, _1, _2)
        base1 = re.sub(r'[_-]?r(ep)?[_-]?[0-9]+', '', stem1)
        base2 = re.sub(r'[_-]?r(ep)?[_-]?[0-9]+', '', stem2)
        
        # Check if bases are similar enough (at least 70% match)
        common_length = sum(1 for a, b in zip(base1, base2) if a == b)
        min_length = min(len(base1), len(base2))
        
        return min_length > 0 and common_length / min_length >= 0.7
    
    @staticmethod
    def _extract_group_id(stem: str) -> str:
        """Extract group identifier (remove rep/replicate numbers)."""
        cleaned = re.sub(r'[_-]?r(ep)?[_-]?[0-9]+', '', stem)
        return cleaned or stem


# Curated palette for BigWig/coverage tracks — visually distinct, colorblind-friendly
_BIGWIG_PALETTE = [
    "#4ECDC4",  # Teal
    "#E8A4D0",  # Rose
    "#95E1D3",  # Mint
    "#F38181",  # Coral
    "#A8D8EA",  # Sky blue
    "#AA96DA",  # Lavender
    "#FCBAD3",  # Pink
    "#FFFFD2",  # Light yellow
]

# Fixed color for annotation-style tracks (narrowpeak, bed)
_ANNOTATION_COLOR = "#FF6B6B"


def _palette_color_for_group(group_name: str) -> str:
    """Deterministically map a group name to a palette color."""
    idx = hash(group_name) % len(_BIGWIG_PALETTE)
    return _BIGWIG_PALETTE[idx]


class InferenceResult:
    """Result of inference operations with confidence and explanations."""

    def __init__(self):
        self.track_type: TemplateTrackType = TemplateTrackType.UNKNOWN
        self.type_confidence: float = 0.0
        self.title: str = "Unknown"
        self.title_inferred: bool = False
        self.group: Optional[str] = None
        self.group_confidence: float = 0.0
        self.suggested_color: Optional[str] = None  # Suggested default color
        self.notes: list[str] = []  # Reasoning notes
        self.issues: list[str] = []  # Warnings or ambiguities

        # Seqnado-specific
        self.is_seqnado: bool = False
        self.seqnado_sample: Optional[str] = None
        self.seqnado_antibody: Optional[str] = None

    def overall_confidence(self) -> float:
        """Average confidence across all inferences."""
        confidences = [self.type_confidence]
        if self.group:
            confidences.append(self.group_confidence)
        return sum(confidences) / len(confidences) if confidences else 0.0


def infer_track(path: str, known_group: Optional[str] = None) -> InferenceResult:
    """
    Run full inference on a single track.
    
    Args:
        path: File path or URL
        known_group: Pre-determined group (overrides heuristic grouping)
    
    Returns:
        InferenceResult with type, title, group, and confidence
    """
    result = InferenceResult()
    
    # Check for seqnado pattern
    seqnado = SeqnadoPattern.parse(Path(path).name)
    if seqnado:
        result.is_seqnado = True
        result.seqnado_sample, result.seqnado_antibody = seqnado
        result.notes.append(
            f"Seqnado pattern detected: sample='{result.seqnado_sample}', "
            f"antibody='{result.seqnado_antibody}'"
        )
    
    # Classify track type
    track_type, type_conf = TrackClassifier.classify(path)
    result.track_type = track_type
    result.type_confidence = type_conf
    
    if type_conf < 0.5:
        result.issues.append(f"Could not confidently determine track type for {path}")
        result.notes.append(f"Using default type: {track_type.value}")
    else:
        result.notes.append(f"Detected track type: {track_type.value} ({type_conf:.0%} confidence)")
    
    # Infer title
    title, was_inferred = TitleInference.infer(path, track_type)
    result.title = title
    result.title_inferred = was_inferred
    result.notes.append(f"Title: {title}" + (" (inferred)" if was_inferred else " (explicit)"))
    
    # Handle grouping
    if known_group:
        result.group = known_group
        result.group_confidence = 1.0
        result.notes.append(f"Assigned to group: {known_group}")

    # Assign suggested color
    if result.track_type in (TemplateTrackType.BIGWIG, TemplateTrackType.BEDGRAPH):
        # Use antibody name for seqnado files (consistent color per antibody across samples)
        color_key = result.seqnado_antibody or result.group or path
        result.suggested_color = _palette_color_for_group(color_key)
    elif result.track_type in (TemplateTrackType.NARROWPEAK, TemplateTrackType.BED, TemplateTrackType.ANNOTATION):
        result.suggested_color = _ANNOTATION_COLOR
    # UNKNOWN, GENE, LINKS, OVERLAY: no color suggestion

    return result
