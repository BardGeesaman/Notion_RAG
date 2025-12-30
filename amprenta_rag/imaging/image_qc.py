"""Image quality control pipeline for microscopy images."""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import ndimage
from scipy.ndimage import gaussian_filter, median_filter
from skimage import filters, measure, morphology


@dataclass
class FocusResult:
    """Result of focus quality assessment."""
    score: float  # Higher is better (0-1 normalized)
    algorithm: str  # "laplacian", "brenner", "variance"
    is_focused: bool  # Above threshold
    raw_value: float  # Un-normalized metric value


@dataclass
class SaturationResult:
    """Result of saturation detection."""
    saturated_percent: float  # 0-100
    saturated_pixel_count: int
    total_pixels: int
    is_saturated: bool  # Above threshold (typically >1%)
    hot_spots: List[Tuple[int, int]] = field(default_factory=list)  # (x, y) of saturated regions


@dataclass
class UniformityResult:
    """Result of illumination uniformity assessment."""
    uniformity_score: float  # 0-1, higher is more uniform
    vignetting_detected: bool
    gradient_direction: Optional[str] = None  # "left-right", "top-bottom", etc.
    corner_ratios: Dict[str, float] = field(default_factory=dict)  # {"TL": 0.8, "TR": 0.9, ...}


@dataclass
class ArtifactResult:
    """Result of artifact detection."""
    artifact_count: int
    artifact_percent: float  # Area covered
    artifact_types: List[str] = field(default_factory=list)  # "dust", "debris", "scratch"
    artifact_locations: List[Tuple[int, int, int, int]] = field(default_factory=list)  # bounding boxes


@dataclass
class ImageQCResult:
    """Complete QC result for a single image."""
    image_path: str
    focus: FocusResult
    saturation: SaturationResult
    uniformity: UniformityResult
    artifacts: Optional[ArtifactResult] = None
    overall_score: float = 0.0  # 0-100 composite
    passed_qc: bool = False
    issues: List[str] = field(default_factory=list)
    timestamp: datetime = field(default_factory=datetime.now)


@dataclass
class PlateQCReport:
    """QC report for an entire plate."""
    plate_id: str
    total_images: int
    passed_count: int
    failed_count: int
    average_focus_score: float
    focus_heatmap: Dict[str, float] = field(default_factory=dict)  # {well_position: score}
    saturation_alerts: List[str] = field(default_factory=list)  # Well positions with saturation
    uniformity_issues: List[str] = field(default_factory=list)  # Well positions with vignetting
    recommendations: List[str] = field(default_factory=list)


def calculate_focus_score(
    image: np.ndarray, 
    algorithm: str = "laplacian",
    normalize: bool = True
) -> FocusResult:
    """
    Calculate image focus quality.
    
    Args:
        image: Input image (grayscale or will be converted)
        algorithm: Focus metric algorithm ("laplacian", "brenner", "variance")
        normalize: Whether to normalize score to 0-1 range
        
    Returns:
        FocusResult with score and metadata
    """
    # Convert to grayscale if needed
    if len(image.shape) == 3:
        # Convert RGB to grayscale using standard weights
        image = 0.299 * image[:, :, 0] + 0.587 * image[:, :, 1] + 0.114 * image[:, :, 2]
    
    # Ensure float type for calculations
    image = image.astype(np.float64)
    
    if algorithm == "laplacian":
        # Variance of Laplacian - widely used and reliable
        laplacian = ndimage.laplace(image)
        raw_value = float(np.var(laplacian))
        
        # Normalize by image intensity range for robustness
        if normalize and image.max() > image.min():
            intensity_range = image.max() - image.min()
            normalized_score = min(1.0, raw_value / (intensity_range ** 2 * 100))
        else:
            normalized_score = min(1.0, raw_value / 10000.0)  # Empirical scaling
            
    elif algorithm == "brenner":
        # Brenner gradient - more robust to noise
        # Sum of squared differences between pixels 2 positions apart
        diff_x = np.diff(image, n=2, axis=1)
        diff_y = np.diff(image, n=2, axis=0)
        raw_value = float(np.sum(diff_x**2) + np.sum(diff_y**2))
        
        if normalize:
            # Normalize by image size and intensity
            pixel_count = image.size
            intensity_range = image.max() - image.min() if image.max() > image.min() else 1.0
            normalized_score = min(1.0, raw_value / (pixel_count * intensity_range**2))
        else:
            normalized_score = raw_value
            
    elif algorithm == "variance":
        # Normalized variance - simple but effective
        raw_value = float(np.var(image))
        
        if normalize:
            # Normalize by intensity range
            intensity_range = image.max() - image.min() if image.max() > image.min() else 1.0
            normalized_score = min(1.0, raw_value / (intensity_range**2))
        else:
            normalized_score = raw_value
            
    else:
        raise ValueError(f"Unknown focus algorithm: {algorithm}")
    
    # Determine if focused (threshold depends on algorithm)
    focus_thresholds = {
        "laplacian": 0.001,  # Lower threshold for normalized scores
        "brenner": 0.001,
        "variance": 0.001
    }
    
    score = normalized_score if normalize else raw_value
    is_focused = score > focus_thresholds.get(algorithm, 0.1)
    
    return FocusResult(
        score=score,
        algorithm=algorithm,
        is_focused=is_focused,
        raw_value=raw_value
    )


def detect_saturation(
    image: np.ndarray,
    threshold: float = 0.99,  # 99% of max value
    min_cluster_size: int = 10
) -> SaturationResult:
    """
    Detect saturated pixels and hot spots.
    
    Args:
        image: Input image
        threshold: Saturation threshold as fraction of max value
        min_cluster_size: Minimum pixels for hot spot detection
        
    Returns:
        SaturationResult with saturation statistics
    """
    # Determine max value based on data type
    if image.dtype == np.uint8:
        max_value = 255
    elif image.dtype == np.uint16:
        max_value = 65535
    else:
        max_value = float(np.max(image))
    
    # Find saturated pixels
    saturation_level = threshold * max_value
    saturated_mask = image >= saturation_level
    
    saturated_pixel_count = int(np.sum(saturated_mask))
    total_pixels = int(image.size)
    saturated_percent = (saturated_pixel_count / total_pixels) * 100.0
    
    # Find hot spots (clusters of saturated pixels)
    hot_spots = []
    if saturated_pixel_count > 0:
        # Label connected components of saturated pixels
        labeled_saturated = measure.label(saturated_mask)
        regions = measure.regionprops(labeled_saturated)
        
        for region in regions:
            if region.area >= min_cluster_size:
                # Use centroid as hot spot location
                y, x = region.centroid
                hot_spots.append((int(x), int(y)))
    
    # Typically consider >1% saturation as problematic
    is_saturated = saturated_percent > 1.0
    
    return SaturationResult(
        saturated_percent=saturated_percent,
        saturated_pixel_count=saturated_pixel_count,
        total_pixels=total_pixels,
        is_saturated=is_saturated,
        hot_spots=hot_spots
    )


def assess_illumination_uniformity(
    image: np.ndarray,
    grid_size: int = 8  # 8x8 grid analysis
) -> UniformityResult:
    """
    Assess illumination uniformity across image.
    
    Args:
        image: Input image
        grid_size: Grid size for uniformity analysis
        
    Returns:
        UniformityResult with uniformity metrics
    """
    # Convert to grayscale if needed
    if len(image.shape) == 3:
        image = 0.299 * image[:, :, 0] + 0.587 * image[:, :, 1] + 0.114 * image[:, :, 2]
    
    height, width = image.shape
    
    # Create grid and calculate mean intensity in each cell
    grid_means = []
    grid_positions = []
    
    cell_height = height // grid_size
    cell_width = width // grid_size
    
    for i in range(grid_size):
        for j in range(grid_size):
            y_start = i * cell_height
            y_end = min((i + 1) * cell_height, height)
            x_start = j * cell_width
            x_end = min((j + 1) * cell_width, width)
            
            cell = image[y_start:y_end, x_start:x_end]
            if cell.size > 0:
                grid_means.append(np.mean(cell))
                grid_positions.append((i, j))
    
    if not grid_means:
        return UniformityResult(
            uniformity_score=0.0,
            vignetting_detected=True,
            gradient_direction=None,
            corner_ratios={}
        )
    
    grid_means = np.array(grid_means)
    
    # Calculate uniformity score (1 - coefficient of variation)
    if np.mean(grid_means) > 0:
        cv = np.std(grid_means) / np.mean(grid_means)
        uniformity_score = max(0.0, 1.0 - cv)
    else:
        uniformity_score = 0.0
    
    # Analyze corner regions vs center for vignetting
    center_i, center_j = grid_size // 2, grid_size // 2
    grid_array = np.array(grid_means).reshape(grid_size, grid_size)
    
    center_intensity = grid_array[center_i, center_j] if grid_size > 2 else np.mean(grid_means)
    
    corner_ratios = {}
    if grid_size >= 3:
        corners = {
            "TL": grid_array[0, 0],           # Top-left
            "TR": grid_array[0, -1],          # Top-right
            "BL": grid_array[-1, 0],          # Bottom-left
            "BR": grid_array[-1, -1]          # Bottom-right
        }
        
        for corner, intensity in corners.items():
            if center_intensity > 0:
                corner_ratios[corner] = float(intensity / center_intensity)
            else:
                corner_ratios[corner] = 1.0
    
    # Detect vignetting (corners significantly darker than center)
    vignetting_threshold = 0.8
    vignetting_detected = any(ratio < vignetting_threshold for ratio in corner_ratios.values())
    
    # Analyze gradient direction
    gradient_direction = None
    if grid_size >= 3:
        left_mean = np.mean(grid_array[:, 0])
        right_mean = np.mean(grid_array[:, -1])
        top_mean = np.mean(grid_array[0, :])
        bottom_mean = np.mean(grid_array[-1, :])
        
        horizontal_diff = abs(left_mean - right_mean)
        vertical_diff = abs(top_mean - bottom_mean)
        
        # Significant gradient if difference > 10% of mean intensity
        mean_intensity = np.mean(grid_means)
        threshold = 0.1 * mean_intensity
        
        if horizontal_diff > threshold and horizontal_diff > vertical_diff:
            gradient_direction = "left-right" if left_mean > right_mean else "right-left"
        elif vertical_diff > threshold:
            gradient_direction = "top-bottom" if top_mean > bottom_mean else "bottom-top"
    
    return UniformityResult(
        uniformity_score=uniformity_score,
        vignetting_detected=vignetting_detected,
        gradient_direction=gradient_direction,
        corner_ratios=corner_ratios
    )


def detect_artifacts(
    image: np.ndarray,
    sensitivity: float = 0.5
) -> ArtifactResult:
    """
    Detect image artifacts (dust, debris, scratches).
    
    Args:
        image: Input image
        sensitivity: Detection sensitivity (0.0-1.0)
        
    Returns:
        ArtifactResult with artifact information
    """
    # Convert to grayscale if needed
    if len(image.shape) == 3:
        image = 0.299 * image[:, :, 0] + 0.587 * image[:, :, 1] + 0.114 * image[:, :, 2]
    
    # Smooth image to get background estimation
    background = gaussian_filter(image, sigma=10)
    
    # Calculate difference from background
    diff = np.abs(image - background)
    
    # Adaptive threshold based on local statistics
    local_median = median_filter(diff, size=15)
    local_std = ndimage.uniform_filter(diff**2, size=15) - ndimage.uniform_filter(diff, size=15)**2
    local_std = np.sqrt(np.maximum(local_std, 0))
    
    # Threshold for artifact detection
    threshold_multiplier = 2.0 + (1.0 - sensitivity) * 3.0  # 2-5x std depending on sensitivity
    threshold = local_median + threshold_multiplier * local_std
    
    # Find potential artifacts
    artifact_mask = diff > threshold
    
    # Remove small noise with morphological operations
    min_artifact_size = max(10, int(image.size * 0.0001))  # At least 0.01% of image
    artifact_mask = morphology.remove_small_objects(artifact_mask, min_size=min_artifact_size)
    
    # Label connected components
    labeled_artifacts = measure.label(artifact_mask)
    regions = measure.regionprops(labeled_artifacts)
    
    artifact_count = len(regions)
    artifact_percent = (np.sum(artifact_mask) / image.size) * 100.0
    
    # Classify artifacts and get bounding boxes
    artifact_types = []
    artifact_locations = []
    
    for region in regions:
        # Get bounding box
        min_row, min_col, max_row, max_col = region.bbox
        artifact_locations.append((min_col, min_row, max_col, max_row))
        
        # Simple classification based on shape
        aspect_ratio = region.major_axis_length / max(region.minor_axis_length, 1e-6)
        eccentricity = region.eccentricity
        
        if aspect_ratio > 5.0:
            artifact_types.append("scratch")
        elif eccentricity < 0.5 and region.area < image.size * 0.01:
            artifact_types.append("dust")
        else:
            artifact_types.append("debris")
    
    return ArtifactResult(
        artifact_count=artifact_count,
        artifact_percent=artifact_percent,
        artifact_types=artifact_types,
        artifact_locations=artifact_locations
    )


def calculate_signal_noise_ratio(image: np.ndarray) -> float:
    """
    Calculate SNR using mean/std of foreground vs background.
    
    Args:
        image: Input image
        
    Returns:
        Signal-to-noise ratio
    """
    # Convert to grayscale if needed
    if len(image.shape) == 3:
        image = 0.299 * image[:, :, 0] + 0.587 * image[:, :, 1] + 0.114 * image[:, :, 2]
    
    # Use Otsu's threshold to separate foreground and background
    try:
        threshold = filters.threshold_otsu(image)
        foreground_mask = image > threshold
        
        if np.sum(foreground_mask) == 0 or np.sum(~foreground_mask) == 0:
            # Fallback: use image statistics
            signal = np.mean(image)
            noise = np.std(image)
        else:
            # Calculate signal as mean of foreground
            signal = np.mean(image[foreground_mask])
            # Calculate noise as std of background
            noise = np.std(image[~foreground_mask])
        
        if noise > 0:
            snr = signal / noise
        else:
            snr = float('inf') if signal > 0 else 0.0
            
    except Exception:
        # Fallback calculation
        signal = np.mean(image)
        noise = np.std(image)
        snr = signal / noise if noise > 0 else 0.0
    
    return float(snr)


def run_image_qc(
    image: np.ndarray,
    image_path: str = "",
    run_artifact_detection: bool = False  # Expensive operation
) -> ImageQCResult:
    """
    Run all QC checks on single image.
    
    Args:
        image: Input image
        image_path: Path to image file (for metadata)
        run_artifact_detection: Whether to run artifact detection (slow)
        
    Returns:
        ImageQCResult with all QC metrics
    """
    # Suppress warnings during QC processing
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        
        # Run individual QC checks
        focus = calculate_focus_score(image, algorithm="laplacian")
        saturation = detect_saturation(image)
        uniformity = assess_illumination_uniformity(image)
        
        artifacts = None
        if run_artifact_detection:
            artifacts = detect_artifacts(image)
        
        # Calculate overall score (0-100)
        focus_weight = 0.4
        saturation_weight = 0.3
        uniformity_weight = 0.3
        
        focus_score = focus.score * 100
        saturation_score = max(0, 100 - saturation.saturated_percent * 10)  # Penalize saturation
        uniformity_score = uniformity.uniformity_score * 100
        
        overall_score = (
            focus_weight * focus_score +
            saturation_weight * saturation_score +
            uniformity_weight * uniformity_score
        )
        
        # Determine if image passes QC
        passed, issues = qc_threshold_check(
            ImageQCResult(
                image_path=image_path,
                focus=focus,
                saturation=saturation,
                uniformity=uniformity,
                artifacts=artifacts,
                overall_score=overall_score,
                passed_qc=False,  # Will be updated
                issues=[],
                timestamp=datetime.now()
            )
        )
    
    return ImageQCResult(
        image_path=image_path,
        focus=focus,
        saturation=saturation,
        uniformity=uniformity,
        artifacts=artifacts,
        overall_score=overall_score,
        passed_qc=passed,
        issues=issues,
        timestamp=datetime.now()
    )


def generate_plate_qc_report(
    images: List[Tuple[str, np.ndarray]],  # (well_position, image)
    plate_id: str
) -> PlateQCReport:
    """
    Generate plate-wide QC report with heatmaps.
    
    Args:
        images: List of (well_position, image) tuples
        plate_id: Plate identifier
        
    Returns:
        PlateQCReport with aggregated metrics
    """
    total_images = len(images)
    passed_count = 0
    focus_scores = []
    focus_heatmap = {}
    saturation_alerts = []
    uniformity_issues = []
    
    for well_position, image in images:
        qc_result = run_image_qc(image, image_path=f"{plate_id}_{well_position}")
        
        if qc_result.passed_qc:
            passed_count += 1
        
        focus_scores.append(qc_result.focus.score)
        focus_heatmap[well_position] = qc_result.focus.score
        
        if qc_result.saturation.is_saturated:
            saturation_alerts.append(well_position)
        
        if qc_result.uniformity.vignetting_detected:
            uniformity_issues.append(well_position)
    
    failed_count = total_images - passed_count
    average_focus_score = float(np.mean(focus_scores)) if focus_scores else 0.0
    
    # Generate recommendations
    recommendations = []
    if failed_count > total_images * 0.1:  # >10% failure rate
        recommendations.append("High failure rate detected - check imaging conditions")
    
    if len(saturation_alerts) > total_images * 0.05:  # >5% saturation
        recommendations.append("Multiple wells show saturation - reduce exposure time")
    
    if len(uniformity_issues) > total_images * 0.2:  # >20% uniformity issues
        recommendations.append("Illumination uniformity issues detected - check light source")
    
    if average_focus_score < 0.3:
        recommendations.append("Overall focus quality is poor - check autofocus settings")
    
    return PlateQCReport(
        plate_id=plate_id,
        total_images=total_images,
        passed_count=passed_count,
        failed_count=failed_count,
        average_focus_score=average_focus_score,
        focus_heatmap=focus_heatmap,
        saturation_alerts=saturation_alerts,
        uniformity_issues=uniformity_issues,
        recommendations=recommendations
    )


def qc_threshold_check(
    result: ImageQCResult,
    focus_threshold: float = 0.3,
    saturation_threshold: float = 1.0,
    uniformity_threshold: float = 0.7
) -> Tuple[bool, List[str]]:
    """
    Check if image passes QC thresholds.
    
    Args:
        result: ImageQCResult to check
        focus_threshold: Minimum focus score (0-1)
        saturation_threshold: Maximum saturation percent
        uniformity_threshold: Minimum uniformity score (0-1)
        
    Returns:
        Tuple of (passed, list_of_issues)
    """
    issues = []
    
    # Check focus
    if result.focus.score < focus_threshold:
        issues.append(f"Poor focus: {result.focus.score:.3f} < {focus_threshold}")
    
    # Check saturation
    if result.saturation.saturated_percent > saturation_threshold:
        issues.append(f"High saturation: {result.saturation.saturated_percent:.1f}% > {saturation_threshold}%")
    
    # Check uniformity
    if result.uniformity.uniformity_score < uniformity_threshold:
        issues.append(f"Poor uniformity: {result.uniformity.uniformity_score:.3f} < {uniformity_threshold}")
    
    # Check for vignetting
    if result.uniformity.vignetting_detected:
        issues.append("Vignetting detected")
    
    passed = len(issues) == 0
    return passed, issues
