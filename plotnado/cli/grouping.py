"""
Grouping strategies for organizing tracks into shared scaling/coloring groups.

Supports both predefined strategies (seqnado-aware) and custom regex patterns.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import re

from plotnado.cli.inference import SeqnadoPattern


@dataclass
class GroupingResult:
    """Result of applying a grouping strategy."""
    
    groups: dict[str, list[int]]  # group_name -> list of track indices
    strategy_name: str  # Name of the strategy used
    explanation: str  # Human-readable explanation of the grouping


class GroupingStrategy(ABC):
    """Abstract base class for grouping strategies."""
    
    @abstractmethod
    def apply(self, paths: list[str]) -> Optional[GroupingResult]:
        """
        Apply the strategy to a list of file paths.
        
        Returns:
            GroupingResult if strategy applies, None if not applicable
        """
        pass


class HistoneMarkStrategy(GroupingStrategy):
    """Group by histone mark (H3K*, H4K*, etc.) detected in filenames."""
    
    # Common histone marks pattern: H3K4me3, H3K27ac, H4K20me1, H2AK119ub, etc.
    HISTONE_MARK_PATTERN = re.compile(
        r'(H[0-9]K[0-9]+(?:me|ac|ph|ub)[0-9]*)',
        re.IGNORECASE
    )
    
    def apply(self, paths: list[str]) -> Optional[GroupingResult]:
        """Group files by histone mark if detected in filenames."""
        marks: dict[str, list[int]] = {}
        
        for i, path in enumerate(paths):
            filename = Path(path).name
            match = self.HISTONE_MARK_PATTERN.search(filename)
            if match:
                mark = match.group(1)  # e.g., "H3K4me3"
                if mark not in marks:
                    marks[mark] = []
                marks[mark].append(i)
        
        # Only return if we found marks in all files
        if len(marks) == 0 or sum(len(v) for v in marks.values()) != len(paths):
            return None
        
        # Only group marks that appear multiple times
        multi_mark_groups = {k: v for k, v in marks.items() if len(v) > 1}
        single_mark_groups = {k: v for k, v in marks.items() if len(v) == 1}
        
        # Combine: multi-track groups first, then single tracks
        all_groups = {**multi_mark_groups, **single_mark_groups}
        
        if not all_groups:
            return None
        
        return GroupingResult(
            groups={f"{name.lower()}_group": indices for name, indices in all_groups.items()},
            strategy_name="histone-mark",
            explanation=f"Grouped by histone marks: {', '.join(sorted(all_groups.keys()))}",
        )


class SeqnadoSampleStrategy(GroupingStrategy):
    """Group by sample name (seqnado SAMPLE_ANTIBODY.bw pattern)."""
    
    def apply(self, paths: list[str]) -> Optional[GroupingResult]:
        """Group seqnado files by sample name."""
        seqnado_files = [SeqnadoPattern.parse(Path(p).name) for p in paths]
        all_seqnado = all(s is not None for s in seqnado_files)
        
        if not all_seqnado:
            return None
        
        samples: dict[str, list[int]] = {}
        for i, (sample, antibody) in enumerate(seqnado_files):
            if sample not in samples:
                samples[sample] = []
            samples[sample].append(i)
        
        # Only return groups with multiple members
        groups = {k: v for k, v in samples.items() if len(v) > 1}
        
        if not groups:
            return None
        
        return GroupingResult(
            groups={f"{name}_autoscale": indices for name, indices in groups.items()},
            strategy_name="seqnado-sample",
            explanation="Grouped by sample name (seqnado pattern: SAMPLE_ANTIBODY.bw)",
        )


class SeqnadoAntibodyStrategy(GroupingStrategy):
    """Group by antibody name (seqnado SAMPLE_ANTIBODY.bw pattern)."""
    
    def apply(self, paths: list[str]) -> Optional[GroupingResult]:
        """Group seqnado files by antibody name."""
        seqnado_files = [SeqnadoPattern.parse(Path(p).name) for p in paths]
        all_seqnado = all(s is not None for s in seqnado_files)
        
        if not all_seqnado:
            return None
        
        antibodies: dict[str, list[int]] = {}
        for i, (sample, antibody) in enumerate(seqnado_files):
            if antibody not in antibodies:
                antibodies[antibody] = []
            antibodies[antibody].append(i)
        
        # Only return groups with multiple members
        groups = {k: v for k, v in antibodies.items() if len(v) > 1}
        
        if not groups:
            return None
        
        return GroupingResult(
            groups={f"{name}_autoscale": indices for name, indices in groups.items()},
            strategy_name="seqnado-antibody",
            explanation="Grouped by antibody name (seqnado pattern: SAMPLE_ANTIBODY.bw)",
        )


class RegexGroupingStrategy(GroupingStrategy):
    """Group files by regex pattern matching."""
    
    def __init__(self, pattern: str):
        """
        Initialize with a regex pattern.
        
        The pattern should contain a capturing group for the group name.
        Examples:
            r'([^_]+)_rep[0-9]+'  # Matches: control_rep1, control_rep2 -> groups as "control"
            r'(.*?)_[0-9]{8}$'    # Matches: sample_20260312 -> group as "sample"
            r'([ACGT]+)\\.bw'      # Matches: ACGTNN.bw -> group as "ACGTNN"
        """
        self.pattern = re.compile(pattern, re.IGNORECASE)
    
    def apply(self, paths: list[str]) -> Optional[GroupingResult]:
        """Group files by regex pattern."""
        # Extract filenames
        filenames = [Path(p).name for p in paths]
        
        # Try to extract group names from filenames
        groups: dict[str, list[int]] = {}
        matched_any = False
        
        for i, filename in enumerate(filenames):
            match = self.pattern.search(filename)
            if match:
                matched_any = True
                if match.groups():
                    # Use the first capturing group as the group name
                    group_name = match.group(1)
                    if group_name not in groups:
                        groups[group_name] = []
                    groups[group_name].append(i)
                else:
                    # If no capturing group, use the entire match
                    group_name = match.group(0)
                    if group_name not in groups:
                        groups[group_name] = []
                    groups[group_name].append(i)
        
        if not matched_any:
            return None
        
        # Only return groups with multiple members
        groups = {k: v for k, v in groups.items() if len(v) > 1}
        
        if not groups:
            return None
        
        return GroupingResult(
            groups={f"{name}_group": indices for name, indices in groups.items()},
            strategy_name=f"regex: {self.pattern.pattern}",
            explanation=f"Grouped by regex pattern: {self.pattern.pattern}",
        )


class PredefinedGroupingStrategies:
    """Registry of predefined grouping strategies."""
    
    STRATEGIES = {
        "sample": SeqnadoSampleStrategy(),
        "antibody": SeqnadoAntibodyStrategy(),
    }
    
    @classmethod
    def get(cls, name: str) -> Optional[GroupingStrategy]:
        """Get a predefined strategy by name."""
        return cls.STRATEGIES.get(name)
    
    @classmethod
    def list_names(cls) -> list[str]:
        """List all predefined strategy names."""
        return list(cls.STRATEGIES.keys())
    
    @classmethod
    def parse_group_by(cls, group_by_arg: str) -> Optional[GroupingStrategy]:
        """
        Parse a --group-by argument.
        
        Can be:
        - A predefined strategy name (e.g., "sample", "antibody")
        - A regex pattern (detected by presence of regex metacharacters)
        
        Returns:
            GroupingStrategy instance or None if invalid
        """
        if not group_by_arg:
            return None
        
        # Try predefined strategies first
        strategy = cls.get(group_by_arg)
        if strategy:
            return strategy
        
        # Try as regex pattern
        try:
            return RegexGroupingStrategy(group_by_arg)
        except re.error as e:
            raise ValueError(f"Invalid regex pattern: {e}")


def apply_grouping_strategy(
    paths: list[str],
    strategy: GroupingStrategy,
) -> Optional[GroupingResult]:
    """
    Apply a grouping strategy to a list of file paths.
    
    Returns:
        GroupingResult if applicable, None if strategy doesn't apply
    """
    return strategy.apply(paths)


def detect_and_apply_grouping(paths: list[str]) -> Optional[GroupingResult]:
    """
    Automatically detect and apply the best grouping strategy.
    
    Tries strategies in order of preference (most specific first):
    1. Seqnado sample grouping (SAMPLE_ANTIBODY.bw)
    2. Seqnado antibody grouping
    3. Histone mark grouping (H3K*, H4K*, etc.)
    """
    strategies_to_try = [
        SeqnadoSampleStrategy(),
        SeqnadoAntibodyStrategy(),
        HistoneMarkStrategy(),
    ]
    
    for strategy in strategies_to_try:
        result = strategy.apply(paths)
        if result:
            return result
    
    return None
