from typing import List, Optional

from pydantic import BaseModel


class OmicsDatasetIngestRequest(BaseModel):
    omics_type: str
    dataset_path: str
    metadata: Optional[dict] = None
    extra_args: Optional[dict] = None


class OmicsFileMeta(BaseModel):
    file_name: str
    file_format: str
    file_size: Optional[int] = None
    checksum: Optional[str] = None


class OmicsFeature(BaseModel):
    feature_id: str
    attributes: dict


class ParsedOmicsDataset(BaseModel):
    dataset_id: str
    features: List[OmicsFeature]
    metadata: dict
    files: Optional[List[OmicsFileMeta]] = None
