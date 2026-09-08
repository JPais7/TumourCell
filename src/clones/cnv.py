"""CNV evidence container."""
from dataclasses import dataclass
@dataclass(frozen=True)
class CNVProfile:
    sample_id: str; segments: tuple[tuple[str,int,int,float], ...]; quality: str
    @property
    def informative(self): return self.quality == "PASS" and bool(self.segments)
