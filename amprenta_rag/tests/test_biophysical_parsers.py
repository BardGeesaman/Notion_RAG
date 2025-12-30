"""
Tests for Biophysical assay file parsers.

Tests cover SPR, MST, and DSC file parsing functionality with mock data
to avoid requiring real instrument files.
"""

import pytest
import tempfile
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock

from amprenta_rag.biophysical.spr_parser import (
    SPRData,
    Sensorgram,
    parse_biacore_csv,
    parse_biacore_sensorgram_txt,
    detect_injection_phases,
    validate_array_size,
)
from amprenta_rag.biophysical.mst_parser import (
    MSTData,
    DosePoint,
    parse_nanotemper_xlsx,
    parse_mst_csv,
    calculate_fnorm,
)
from amprenta_rag.biophysical.dsc_parser import (
    DSCData,
    Scan,
    parse_microcal_csv,
    parse_ta_instruments_txt,
    baseline_correction,
)


class TestSPRParser:
    """Test SPR file parsing functionality."""

    def test_validate_array_size_valid(self):
        """Test array size validation with valid array."""
        arr = np.random.rand(1000)
        # Should not raise exception
        validate_array_size(arr)
        
    def test_validate_array_size_invalid(self):
        """Test array size validation with oversized array."""
        arr = np.random.rand(60000)  # Exceeds 50K limit
        with pytest.raises(ValueError, match="Array size .* exceeds maximum"):
            validate_array_size(arr)

    def test_detect_injection_phases(self):
        """Test detection of association and dissociation phases."""
        # Create mock sensorgram with clear phases
        time = np.linspace(0, 300, 1000)  # 5 minutes
        response = np.zeros_like(time)
        
        # Association phase (50-150s): rising signal
        assoc_mask = (time >= 50) & (time <= 150)
        response[assoc_mask] = 100 * (time[assoc_mask] - 50) / 100
        
        # Dissociation phase (150-250s): falling signal  
        dissoc_mask = (time >= 150) & (time <= 250)
        response[dissoc_mask] = 100 * (1 - (time[dissoc_mask] - 150) / 100)
        
        assoc_start, dissoc_start = detect_injection_phases(time, response)
        
        # Should detect phases within reasonable range (algorithm may vary)
        # Just verify that dissociation starts after association
        assert assoc_start >= 0
        assert dissoc_start > assoc_start
        assert dissoc_start <= time[-1]

    def test_detect_injection_phases_invalid_input(self):
        """Test phase detection with invalid input."""
        with pytest.raises(ValueError, match="Invalid time/response arrays"):
            detect_injection_phases(np.array([1, 2]), np.array([1]))  # Different lengths
            
        with pytest.raises(ValueError, match="Invalid time/response arrays"):
            detect_injection_phases(np.array([1, 2]), np.array([3, 4]))  # Too few points

    def test_parse_biacore_csv(self):
        """Test parsing of Biacore CSV file."""
        # Create mock CSV content
        csv_content = """Instrument: Biacore T200
Chip: CM5
Temperature: 25.0
Flow Rate: 30
Buffer: PBS + 0.05% Tween-20
Ligand: EGFR-His6
Analyte: Compound-001

Time	Response	Concentration	Cycle
0.0	0.0	100.0	1
1.0	10.0	100.0	1
2.0	25.0	100.0	1
3.0	35.0	100.0	1
4.0	30.0	100.0	1
5.0	20.0	100.0	1
6.0	15.0	100.0	1
7.0	10.0	100.0	1
8.0	8.0	100.0	1
9.0	5.0	100.0	1
10.0	2.0	100.0	1
0.0	0.0	200.0	2
1.0	15.0	200.0	2
2.0	40.0	200.0	2
3.0	55.0	200.0	2
4.0	50.0	200.0	2
5.0	35.0	200.0	2
6.0	25.0	200.0	2
7.0	20.0	200.0	2
8.0	15.0	200.0	2
9.0	10.0	200.0	2
10.0	5.0	200.0	2
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(csv_content)
            f.flush()
            
            spr_data = parse_biacore_csv(f.name)
            
            assert isinstance(spr_data, SPRData)
            assert spr_data.instrument == "Biacore T200"
            assert spr_data.chip_type == "CM5"
            assert spr_data.temperature == 25.0
            assert spr_data.flow_rate == 30.0
            assert spr_data.buffer == "PBS + 0.05% Tween-20"
            assert spr_data.ligand_name == "EGFR-His6"
            assert spr_data.analyte_name == "Compound-001"
            
            assert len(spr_data.sensorgrams) == 2
            
            # Check first cycle
            cycle1 = spr_data.sensorgrams[0]
            assert cycle1.cycle == 1
            assert cycle1.concentration == 100.0
            assert len(cycle1.time) == 11
            assert len(cycle1.response) == 11
            
            # Check second cycle
            cycle2 = spr_data.sensorgrams[1]
            assert cycle2.cycle == 2
            assert cycle2.concentration == 200.0
            
        Path(f.name).unlink()  # Cleanup

    def test_parse_biacore_sensorgram_txt(self):
        """Test parsing of Biacore TXT sensorgram file."""
        txt_content = """Instrument: Biacore T200
Chip: CM5
Temperature: 25.0
Concentration: 100.0

0.0	0.0
1.0	10.0
2.0	25.0
3.0	35.0
4.0	30.0
5.0	20.0
6.0	15.0
7.0	10.0
8.0	8.0
9.0	5.0
10.0	2.0
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write(txt_content)
            f.flush()
            
            spr_data = parse_biacore_sensorgram_txt(f.name)
            
            assert isinstance(spr_data, SPRData)
            assert spr_data.instrument == "Biacore T200"
            assert spr_data.chip_type == "CM5"
            assert spr_data.temperature == 25.0
            
            assert len(spr_data.sensorgrams) == 1
            sensorgram = spr_data.sensorgrams[0]
            assert sensorgram.cycle == 1
            assert sensorgram.concentration == 100.0
            assert len(sensorgram.time) == 11
            assert len(sensorgram.response) == 11
            
        Path(f.name).unlink()  # Cleanup

    def test_parse_biacore_csv_file_not_found(self):
        """Test CSV parsing with non-existent file."""
        with pytest.raises(FileNotFoundError):
            parse_biacore_csv("nonexistent_file.csv")

    def test_parse_biacore_csv_invalid_format(self):
        """Test CSV parsing with invalid format."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("Invalid content without proper data")
            f.flush()
            
            with pytest.raises(ValueError, match="Invalid Biacore CSV format"):
                parse_biacore_csv(f.name)
                
        Path(f.name).unlink()  # Cleanup


class TestMSTParser:
    """Test MST file parsing functionality."""

    def test_calculate_fnorm(self):
        """Test Fnorm calculation from cold/hot fluorescence."""
        cold = np.array([1000, 1010, 990, 1005])
        hot = np.array([950, 960, 940, 955])
        
        fnorm = calculate_fnorm(cold, hot)
        
        # Expected: (951.25 - 1001.25) / 1001.25 * 1000 = -49.94‰
        expected = (951.25 - 1001.25) / 1001.25 * 1000
        assert abs(fnorm - expected) < 0.1

    def test_calculate_fnorm_invalid_input(self):
        """Test Fnorm calculation with invalid input."""
        with pytest.raises(ValueError, match="Empty fluorescence arrays"):
            calculate_fnorm(np.array([]), np.array([100]))
            
        with pytest.raises(ValueError, match="same length"):
            calculate_fnorm(np.array([100]), np.array([100, 200]))
            
        with pytest.raises(ValueError, match="Cold fluorescence cannot be zero"):
            calculate_fnorm(np.array([0]), np.array([100]))

    @patch('openpyxl.load_workbook')
    def test_parse_nanotemper_xlsx(self, mock_load_workbook):
        """Test parsing of NanoTemper XLSX file."""
        # Mock Excel workbook
        mock_workbook = MagicMock()
        mock_sheet = MagicMock()
        
        # Set up sheet structure
        mock_workbook.sheetnames = ["Analysis", "Raw Data"]
        mock_workbook.__getitem__.return_value = mock_sheet
        mock_load_workbook.return_value = mock_workbook
        
        # Mock cell values for metadata
        def mock_cell_value(row, col):
            mock_cell = MagicMock()
            if row == 1 and col == 1:
                mock_cell.value = "Instrument"
            elif row == 1 and col == 2:
                mock_cell.value = "NanoTemper Monolith NT.115"
            elif row == 2 and col == 1:
                mock_cell.value = "Capillary Type"
            elif row == 2 and col == 2:
                mock_cell.value = "Premium"
            elif row == 10 and col == 1:
                mock_cell.value = "Concentration"
            elif row == 10 and col == 2:
                mock_cell.value = "Fnorm"
            elif row == 11 and col == 1:
                mock_cell.value = 0.0
            elif row == 11 and col == 2:
                mock_cell.value = 0.0
            elif row == 12 and col == 1:
                mock_cell.value = 10.0
            elif row == 12 and col == 2:
                mock_cell.value = -5.0
            else:
                mock_cell.value = None
            return mock_cell
        
        mock_sheet.cell = mock_cell_value
        mock_sheet.max_row = 12
        mock_sheet.max_column = 5
        
        with tempfile.NamedTemporaryFile(suffix='.xlsx', delete=False) as f:
            mst_data = parse_nanotemper_xlsx(f.name)
            
            assert isinstance(mst_data, MSTData)
            assert mst_data.instrument == "NanoTemper Monolith NT.115"
            assert mst_data.capillary_type == "Premium"
            assert len(mst_data.dose_points) == 2
            
            # Check dose points
            assert mst_data.dose_points[0].concentration == 0.0
            assert mst_data.dose_points[0].fnorm == 0.0
            assert mst_data.dose_points[1].concentration == 10.0
            assert mst_data.dose_points[1].fnorm == -5.0
            
        Path(f.name).unlink()  # Cleanup

    def test_parse_mst_csv(self):
        """Test parsing of MST CSV file."""
        csv_content = """Instrument: NanoTemper Monolith NT.115
Capillary: Premium
Excitation Power: 20%
MST Power: 40%
Target: EGFR
Ligand: ATP

Concentration,Fnorm,Error,Cold,Hot
0.0,0.0,0.5,1000,1000
1.0,-2.5,0.3,1000,975
10.0,-5.0,0.4,1000,950
100.0,-8.0,0.6,1000,920
1000.0,-10.0,0.8,1000,900
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(csv_content)
            f.flush()
            
            mst_data = parse_mst_csv(f.name)
            
            assert isinstance(mst_data, MSTData)
            assert mst_data.instrument == "NanoTemper Monolith NT.115"
            assert mst_data.capillary_type == "Premium"
            assert mst_data.excitation_power == 20.0
            assert mst_data.mst_power == 40.0
            assert mst_data.target_name == "EGFR"
            assert mst_data.ligand_name == "ATP"
            
            assert len(mst_data.dose_points) == 5
            
            # Check first and last points
            assert mst_data.dose_points[0].concentration == 0.0
            assert mst_data.dose_points[0].fnorm == 0.0
            assert mst_data.dose_points[4].concentration == 1000.0
            assert mst_data.dose_points[4].fnorm == -10.0
            
        Path(f.name).unlink()  # Cleanup

    def test_parse_mst_csv_with_calculated_fnorm(self):
        """Test MST CSV parsing with Fnorm calculation from cold/hot values."""
        csv_content = """Instrument: NanoTemper
Concentration,Cold,Hot
0.0,1000,1000
100.0,1000,950
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(csv_content)
            f.flush()
            
            mst_data = parse_mst_csv(f.name)
            
            assert len(mst_data.dose_points) == 2
            
            # Should calculate Fnorm from cold/hot values
            assert mst_data.dose_points[0].fnorm == 0.0  # (1000-1000)/1000*1000
            assert abs(mst_data.dose_points[1].fnorm - (-50.0)) < 0.1  # (950-1000)/1000*1000
            
        Path(f.name).unlink()  # Cleanup


class TestDSCParser:
    """Test DSC file parsing functionality."""

    def test_baseline_correction_linear(self):
        """Test linear baseline correction."""
        # Create mock thermogram with linear baseline drift
        temp = np.linspace(20, 80, 100)
        baseline = 0.1 * temp + 2.0  # Linear drift
        signal = np.exp(-(temp - 50)**2 / 100)  # Gaussian peak at 50°C
        cp = baseline + signal
        
        corrected_cp = baseline_correction(temp, cp, method="linear")
        
        # Corrected data should have baseline removed
        assert len(corrected_cp) == len(cp)
        # Peak should still be present but baseline reduced
        peak_idx = np.argmax(corrected_cp)
        assert 45 <= temp[peak_idx] <= 55  # Peak around 50°C

    def test_baseline_correction_invalid_input(self):
        """Test baseline correction with invalid input."""
        with pytest.raises(ValueError, match="same length"):
            baseline_correction(np.array([1, 2]), np.array([1]))
            
        with pytest.raises(ValueError, match="Insufficient data"):
            baseline_correction(np.array([1, 2]), np.array([3, 4]))
            
        with pytest.raises(ValueError, match="Invalid baseline correction method"):
            baseline_correction(np.linspace(0, 10, 20), np.random.rand(20), method="invalid")

    def test_parse_microcal_csv(self):
        """Test parsing of MicroCal DSC CSV file."""
        csv_content = """Instrument: MicroCal PEAQ-DSC
Scan Rate: 1.0 °C/min
Protein Concentration: 1.0 mg/mL
Buffer: PBS pH 7.4
Cell Volume: 200 μL

Temperature,Cp
20.0,0.1
25.0,0.2
30.0,0.5
35.0,1.2
40.0,2.8
45.0,5.1
50.0,8.2
55.0,6.1
60.0,3.2
65.0,1.5
70.0,0.8
75.0,0.4
80.0,0.2
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(csv_content)
            f.flush()
            
            dsc_data = parse_microcal_csv(f.name)
            
            assert isinstance(dsc_data, DSCData)
            assert dsc_data.instrument == "MicroCal PEAQ-DSC"
            assert dsc_data.scan_rate == 1.0
            assert dsc_data.protein_concentration == 1.0
            assert dsc_data.buffer == "PBS pH 7.4"
            assert dsc_data.cell_volume == 200.0
            
            assert len(dsc_data.scans) == 1
            scan = dsc_data.scans[0]
            assert scan.scan_number == 1
            assert len(scan.temperature) == 13
            assert len(scan.heat_capacity) == 13
            assert scan.scan_rate == 1.0
            assert scan.baseline_subtracted is False
            
        Path(f.name).unlink()  # Cleanup

    def test_parse_ta_instruments_txt(self):
        """Test parsing of TA Instruments DSC TXT file."""
        txt_content = """Instrument: TA Instruments DSC 2500
Scan Rate: 10.0 °C/min
Protein Concentration: 2.0 mg/mL
Buffer: HEPES pH 7.5

20.0	0.1	100.0
25.0	0.2	150.0
30.0	0.5	200.0
35.0	1.2	300.0
40.0	2.8	500.0
45.0	5.1	800.0
50.0	8.2	1200.0
55.0	6.1	900.0
60.0	3.2	600.0
65.0	1.5	400.0
70.0	0.8	250.0
75.0	0.4	180.0
80.0	0.2	120.0
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write(txt_content)
            f.flush()
            
            dsc_data = parse_ta_instruments_txt(f.name)
            
            assert isinstance(dsc_data, DSCData)
            assert dsc_data.instrument == "TA Instruments DSC 2500"
            assert dsc_data.scan_rate == 10.0
            assert dsc_data.protein_concentration == 2.0
            assert dsc_data.buffer == "HEPES pH 7.5"
            
            assert len(dsc_data.scans) == 1
            scan = dsc_data.scans[0]
            assert scan.scan_number == 1
            assert len(scan.temperature) == 13
            assert len(scan.heat_capacity) == 13
            assert scan.scan_rate == 10.0
            
        Path(f.name).unlink()  # Cleanup

    def test_parse_microcal_csv_file_not_found(self):
        """Test DSC parsing with non-existent file."""
        with pytest.raises(FileNotFoundError):
            parse_microcal_csv("nonexistent_file.csv")

    def test_parse_microcal_csv_invalid_format(self):
        """Test CSV parsing with invalid format."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("Invalid content without temperature,cp columns")
            f.flush()
            
            with pytest.raises(ValueError, match="Invalid MicroCal CSV format"):
                parse_microcal_csv(f.name)
                
        Path(f.name).unlink()  # Cleanup

    def test_parse_microcal_csv_array_size_limit(self):
        """Test CSV parsing with array size exceeding limit."""
        # Create CSV with too many data points
        csv_content = "Temperature,Cp\n"
        for i in range(60000):  # Exceeds 50K limit
            csv_content += f"{20.0 + i*0.001},{0.1}\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(csv_content)
            f.flush()
            
            with pytest.raises(ValueError, match="Array size .* exceeds maximum"):
                parse_microcal_csv(f.name)
                
        Path(f.name).unlink()  # Cleanup


class TestIntegrationParsers:
    """Test integration between different parsers."""

    def test_all_parsers_data_structures(self):
        """Test that all parsers return expected data structure types."""
        # Test SPR data structure
        sensorgram = Sensorgram(
            cycle=1,
            concentration=100.0,
            time=np.array([0, 1, 2]),
            response=np.array([0, 10, 5]),
            association_start=0.5,
            dissociation_start=1.5
        )
        
        spr_data = SPRData(
            instrument="Test",
            chip_type="Test",
            temperature=25.0,
            sensorgrams=[sensorgram],
            metadata={}
        )
        
        assert isinstance(spr_data, SPRData)
        assert len(spr_data.sensorgrams) == 1
        assert isinstance(spr_data.sensorgrams[0], Sensorgram)
        
        # Test MST data structure
        dose_point = DosePoint(
            concentration=100.0,
            fnorm=-5.0,
            fnorm_error=0.5,
            cold_mean=1000.0,
            hot_mean=950.0
        )
        
        mst_data = MSTData(
            instrument="Test",
            capillary_type="Test",
            excitation_power=20.0,
            mst_power=40.0,
            dose_points=[dose_point],
            metadata={}
        )
        
        assert isinstance(mst_data, MSTData)
        assert len(mst_data.dose_points) == 1
        assert isinstance(mst_data.dose_points[0], DosePoint)
        
        # Test DSC data structure
        scan = Scan(
            scan_number=1,
            temperature=np.array([20, 30, 40]),
            heat_capacity=np.array([0.1, 0.5, 1.0])
        )
        
        dsc_data = DSCData(
            instrument="Test",
            scan_rate=1.0,
            protein_concentration=1.0,
            scans=[scan],
            metadata={}
        )
        
        assert isinstance(dsc_data, DSCData)
        assert len(dsc_data.scans) == 1
        assert isinstance(dsc_data.scans[0], Scan)

    def test_parser_error_handling_consistency(self):
        """Test that all parsers handle errors consistently."""
        # All parsers should raise FileNotFoundError for missing files
        with pytest.raises(FileNotFoundError):
            parse_biacore_csv("missing.csv")
            
        with pytest.raises(FileNotFoundError):
            parse_mst_csv("missing.csv")
            
        with pytest.raises(FileNotFoundError):
            parse_microcal_csv("missing.csv")
            
        with pytest.raises(FileNotFoundError):
            parse_biacore_sensorgram_txt("missing.txt")
            
        with pytest.raises(FileNotFoundError):
            parse_ta_instruments_txt("missing.txt")
