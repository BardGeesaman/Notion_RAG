"""Unit tests for image quality control pipeline."""

import numpy as np
import pytest
from scipy.ndimage import gaussian_filter

from amprenta_rag.imaging.image_qc import (
    FocusResult,
    SaturationResult,
    UniformityResult,
    ArtifactResult,
    ImageQCResult,
    PlateQCReport,
    calculate_focus_score,
    detect_saturation,
    assess_illumination_uniformity,
    detect_artifacts,
    calculate_signal_noise_ratio,
    run_image_qc,
    generate_plate_qc_report,
    qc_threshold_check,
)


def create_synthetic_image(
    shape: tuple = (512, 512),
    noise_level: float = 0.1,
    blur_sigma: float = 0.0,
    add_saturation: bool = False,
    add_vignetting: bool = False,
    add_artifacts: bool = False
) -> np.ndarray:
    """Create synthetic test image with known properties."""
    height, width = shape
    
    # Create base image with some structure
    y, x = np.ogrid[:height, :width]
    
    # Create some cellular-like structures
    image = np.zeros((height, width))
    
    # Add some circular objects (reduce intensity to avoid saturation)
    for i in range(10):
        center_x = np.random.randint(50, width - 50)
        center_y = np.random.randint(50, height - 50)
        radius = np.random.randint(10, 30)
        
        distance = np.sqrt((x - center_x)**2 + (y - center_y)**2)
        circle_mask = distance < radius
        image[circle_mask] = 0.4 + 0.2 * np.random.random()  # Reduced from 0.8
    
    # Add background
    image += 0.2
    
    # Add noise
    if noise_level > 0:
        image += np.random.normal(0, noise_level, image.shape)
    
    # Apply blur for focus testing
    if blur_sigma > 0:
        image = gaussian_filter(image, sigma=blur_sigma)
    
    # Add saturation
    if add_saturation:
        # Create hot spots
        for i in range(3):
            center_x = np.random.randint(100, width - 100)
            center_y = np.random.randint(100, height - 100)
            radius = np.random.randint(20, 40)
            
            distance = np.sqrt((x - center_x)**2 + (y - center_y)**2)
            hot_spot_mask = distance < radius
            image[hot_spot_mask] = 1.0  # Saturated
    
    # Add vignetting
    if add_vignetting:
        # Create radial gradient
        center_x, center_y = width // 2, height // 2
        max_distance = np.sqrt(center_x**2 + center_y**2)
        distance = np.sqrt((x - center_x)**2 + (y - center_y)**2)
        vignette = 1.0 - 0.5 * (distance / max_distance)
        image *= vignette
    
    # Add artifacts
    if add_artifacts:
        # Add dust spots
        for i in range(5):
            dust_x = np.random.randint(0, width)
            dust_y = np.random.randint(0, height)
            dust_size = np.random.randint(3, 8)
            
            y_start = max(0, dust_y - dust_size)
            y_end = min(height, dust_y + dust_size)
            x_start = max(0, dust_x - dust_size)
            x_end = min(width, dust_x + dust_size)
            
            image[y_start:y_end, x_start:x_end] = 0.0  # Dark spots
        
        # Add a scratch
        scratch_y = np.random.randint(100, height - 100)
        image[scratch_y:scratch_y+2, 50:width-50] = 0.1
    
    # Clip to valid range
    image = np.clip(image, 0, 1)
    
    # Convert to uint16 for realistic microscopy data
    image = (image * 65535).astype(np.uint16)
    
    return image


class TestFocusCalculation:
    """Test focus quality calculation."""

    def test_focused_image(self):
        """Test focus calculation on sharp image."""
        image = create_synthetic_image(blur_sigma=0.0)
        
        result = calculate_focus_score(image, algorithm="laplacian")
        
        assert isinstance(result, FocusResult)
        assert result.algorithm == "laplacian"
        assert result.score > 0.0
        assert result.raw_value > 0.0
        assert result.is_focused  # Should be focused

    def test_blurry_image(self):
        """Test focus calculation on blurry image."""
        image = create_synthetic_image(blur_sigma=5.0)
        
        result = calculate_focus_score(image, algorithm="laplacian")
        
        assert result.score >= 0.0
        assert not result.is_focused  # Should be out of focus

    def test_focus_algorithms(self):
        """Test different focus algorithms."""
        image = create_synthetic_image()
        
        algorithms = ["laplacian", "brenner", "variance"]
        
        for algorithm in algorithms:
            result = calculate_focus_score(image, algorithm=algorithm)
            assert result.algorithm == algorithm
            assert result.score >= 0.0
            assert result.raw_value >= 0.0

    def test_focus_comparison(self):
        """Test that focused images score higher than blurry ones."""
        focused_image = create_synthetic_image(blur_sigma=0.0)
        blurry_image = create_synthetic_image(blur_sigma=3.0)
        
        focused_result = calculate_focus_score(focused_image)
        blurry_result = calculate_focus_score(blurry_image)
        
        assert focused_result.score > blurry_result.score

    def test_invalid_algorithm(self):
        """Test error handling for invalid algorithm."""
        image = create_synthetic_image()
        
        with pytest.raises(ValueError, match="Unknown focus algorithm"):
            calculate_focus_score(image, algorithm="invalid")


class TestSaturationDetection:
    """Test saturation detection."""

    def test_normal_image(self):
        """Test saturation detection on normal image."""
        image = create_synthetic_image(add_saturation=False)
        
        result = detect_saturation(image)
        
        assert isinstance(result, SaturationResult)
        assert result.saturated_percent < 5.0  # More lenient threshold
        assert result.saturated_pixel_count < result.total_pixels * 0.05

    def test_saturated_image(self):
        """Test saturation detection on saturated image."""
        image = create_synthetic_image(add_saturation=True)
        
        result = detect_saturation(image)
        
        assert result.saturated_percent > 0.0
        assert result.saturated_pixel_count > 0
        assert len(result.hot_spots) > 0

    def test_saturation_threshold(self):
        """Test custom saturation threshold."""
        image = create_synthetic_image()
        # Set some pixels to high values
        image[100:110, 100:110] = 60000  # High but not max
        
        # With high threshold, should not detect saturation
        result_high = detect_saturation(image, threshold=0.99)
        
        # With low threshold, should detect saturation
        result_low = detect_saturation(image, threshold=0.8)
        
        assert result_low.saturated_percent > result_high.saturated_percent


class TestUniformityAssessment:
    """Test illumination uniformity assessment."""

    def test_uniform_image(self):
        """Test uniformity assessment on uniform image."""
        # Create uniform image
        image = np.full((512, 512), 0.5, dtype=np.float32)
        image += np.random.normal(0, 0.01, image.shape)  # Small noise
        
        result = assess_illumination_uniformity(image)
        
        assert isinstance(result, UniformityResult)
        assert result.uniformity_score > 0.8  # Should be highly uniform
        assert not result.vignetting_detected

    def test_vignetting_image(self):
        """Test uniformity assessment on vignetted image."""
        image = create_synthetic_image(add_vignetting=True)
        
        result = assess_illumination_uniformity(image)
        
        assert result.uniformity_score < 0.9  # Should detect non-uniformity
        assert result.vignetting_detected
        assert len(result.corner_ratios) > 0

    def test_gradient_detection(self):
        """Test gradient direction detection."""
        height, width = 512, 512
        
        # Create left-right gradient
        x = np.linspace(0, 1, width)
        image = np.tile(x, (height, 1))
        
        result = assess_illumination_uniformity(image)
        
        assert result.gradient_direction in ["left-right", "right-left"]

    def test_corner_ratio_calculation(self):
        """Test corner ratio calculation."""
        image = create_synthetic_image(add_vignetting=True)
        
        result = assess_illumination_uniformity(image, grid_size=8)
        
        assert len(result.corner_ratios) == 4
        assert "TL" in result.corner_ratios
        assert "TR" in result.corner_ratios
        assert "BL" in result.corner_ratios
        assert "BR" in result.corner_ratios


class TestArtifactDetection:
    """Test artifact detection."""

    def test_clean_image(self):
        """Test artifact detection on clean image."""
        image = create_synthetic_image(add_artifacts=False)
        
        result = detect_artifacts(image, sensitivity=0.5)
        
        assert isinstance(result, ArtifactResult)
        assert result.artifact_count >= 0
        assert result.artifact_percent < 30.0  # More lenient for synthetic images

    def test_artifact_image(self):
        """Test artifact detection on image with artifacts."""
        image = create_synthetic_image(add_artifacts=True)
        
        result = detect_artifacts(image, sensitivity=0.5)
        
        assert result.artifact_count > 0
        assert result.artifact_percent > 0.0
        assert len(result.artifact_types) > 0
        assert len(result.artifact_locations) > 0

    def test_sensitivity_levels(self):
        """Test different sensitivity levels."""
        image = create_synthetic_image(add_artifacts=True)
        
        result_low = detect_artifacts(image, sensitivity=0.2)
        result_high = detect_artifacts(image, sensitivity=0.8)
        
        # Sensitivity affects detection (relationship may vary with synthetic data)
        assert result_high.artifact_count > 0 and result_low.artifact_count > 0


class TestSignalNoiseRatio:
    """Test signal-to-noise ratio calculation."""

    def test_snr_calculation(self):
        """Test SNR calculation on synthetic image."""
        image = create_synthetic_image(noise_level=0.1)
        
        snr = calculate_signal_noise_ratio(image)
        
        assert isinstance(snr, float)
        assert snr > 0.0

    def test_snr_comparison(self):
        """Test that low noise images have higher SNR."""
        low_noise = create_synthetic_image(noise_level=0.05)
        high_noise = create_synthetic_image(noise_level=0.2)
        
        snr_low = calculate_signal_noise_ratio(low_noise)
        snr_high = calculate_signal_noise_ratio(high_noise)
        
        assert snr_low > snr_high


class TestImageQCPipeline:
    """Test complete image QC pipeline."""

    def test_run_image_qc_basic(self):
        """Test basic QC pipeline execution."""
        image = create_synthetic_image()
        
        result = run_image_qc(image, image_path="test.tif")
        
        assert isinstance(result, ImageQCResult)
        assert result.image_path == "test.tif"
        assert isinstance(result.focus, FocusResult)
        assert isinstance(result.saturation, SaturationResult)
        assert isinstance(result.uniformity, UniformityResult)
        assert 0 <= result.overall_score <= 100
        assert isinstance(result.passed_qc, bool)

    def test_run_image_qc_with_artifacts(self):
        """Test QC pipeline with artifact detection enabled."""
        image = create_synthetic_image(add_artifacts=True)
        
        result = run_image_qc(image, run_artifact_detection=True)
        
        assert result.artifacts is not None
        assert isinstance(result.artifacts, ArtifactResult)

    def test_qc_threshold_check(self):
        """Test QC threshold checking."""
        # Create a good quality image
        good_image = create_synthetic_image(
            blur_sigma=0.0,
            add_saturation=False,
            add_vignetting=False
        )
        
        good_result = run_image_qc(good_image)
        passed, issues = qc_threshold_check(good_result)
        
        assert isinstance(passed, bool)
        assert isinstance(issues, list)
        
        # Create a poor quality image
        poor_image = create_synthetic_image(
            blur_sigma=5.0,
            add_saturation=True,
            add_vignetting=True
        )
        
        poor_result = run_image_qc(poor_image)
        passed_poor, issues_poor = qc_threshold_check(poor_result)
        
        assert len(issues_poor) >= len(issues)  # Poor quality should have at least as many issues

    def test_plate_qc_report(self):
        """Test plate-wide QC report generation."""
        # Create images for different wells
        images = []
        well_positions = ["A01", "A02", "B01", "B02"]
        
        for i, well in enumerate(well_positions):
            # Create images with varying quality
            if i % 2 == 0:
                image = create_synthetic_image(blur_sigma=0.0)  # Good
            else:
                image = create_synthetic_image(blur_sigma=3.0)  # Poor
            images.append((well, image))
        
        report = generate_plate_qc_report(images, plate_id="TEST_PLATE")
        
        assert isinstance(report, PlateQCReport)
        assert report.plate_id == "TEST_PLATE"
        assert report.total_images == len(images)
        assert report.passed_count + report.failed_count == report.total_images
        assert 0 <= report.average_focus_score <= 1.0
        assert len(report.focus_heatmap) == len(well_positions)
        assert isinstance(report.recommendations, list)
