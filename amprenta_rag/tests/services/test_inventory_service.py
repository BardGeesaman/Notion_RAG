"""Unit tests for compound inventory service."""

import pytest
from datetime import date, timedelta
from uuid import uuid4
from unittest.mock import MagicMock, patch

from amprenta_rag.services import inventory as service


class TestBarcodeGeneration:
    """Test barcode generation."""
    
    def test_generate_barcode_format(self):
        """Barcode follows PREFIX-YYYYMMDD-XXXXXX format."""
        barcode = service.generate_barcode("COMP")
        parts = barcode.split("-")
        assert len(parts) == 3
        assert parts[0] == "COMP"
        assert len(parts[1]) == 8  # YYYYMMDD
        assert len(parts[2]) == 6  # hex

    def test_generate_barcode_unique(self):
        """Consecutive barcodes are unique."""
        bc1 = service.generate_barcode()
        bc2 = service.generate_barcode()
        assert bc1 != bc2

    def test_generate_barcode_custom_prefix(self):
        """Custom prefix is used in barcode."""
        barcode = service.generate_barcode("CUSTOM")
        assert barcode.startswith("CUSTOM-")


@patch('amprenta_rag.services.inventory.Compound')
@patch('amprenta_rag.services.inventory.Sample')
class TestCompoundSampleCRUD:
    """Test compound sample CRUD operations."""
    
    def test_create_compound_sample_validates_compound(self, mock_sample, mock_compound):
        """Create sample validates compound exists."""
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = None  # No compound
        
        with pytest.raises(ValueError, match="not found"):
            service.create_compound_sample(
                db=mock_db,
                compound_id=uuid4(),
                quantity=100.0,
            )

    def test_create_sample_validates_barcode_uniqueness(self, mock_sample, mock_compound):
        """P1 TEST: Barcode collision raises ValueError."""
        mock_db = MagicMock()
        
        # Mock existing sample with same barcode
        mock_existing_sample = MagicMock()
        mock_db.query.return_value.filter.return_value.first.side_effect = [
            mock_existing_sample,  # Barcode check finds existing
        ]
        
        with pytest.raises(ValueError, match="already in use"):
            service.create_compound_sample(
                db=mock_db,
                compound_id=uuid4(),
                quantity=100.0,
                barcode="TEST-BARCODE",
            )

    def test_create_sample_success(self, mock_sample, mock_compound):
        """Create sample succeeds with valid inputs."""
        mock_db = MagicMock()
        mock_compound_obj = MagicMock()
        mock_compound_obj.compound_id = "TEST-001"
        
        # Mock successful lookups
        mock_db.query.return_value.filter.return_value.first.side_effect = [
            None,  # No existing barcode
            mock_compound_obj,  # Compound exists
        ]
        
        mock_sample_obj = MagicMock()
        mock_sample_obj.id = uuid4()
        mock_sample_obj.barcode = "TEST-BARCODE"
        mock_sample.return_value = mock_sample_obj
        
        result = service.create_compound_sample(
            db=mock_db,
            compound_id=uuid4(),
            quantity=100.0,
            barcode="TEST-BARCODE",
        )
        
        # Verify sample was created and added
        mock_db.add.assert_called_once()
        mock_db.commit.assert_called_once()
        assert result == mock_sample_obj


class TestRequestWorkflow:
    """Test compound request workflow."""
    
    def test_create_request_requires_target(self):
        """Request without sample_id or compound_id raises ValueError."""
        mock_db = MagicMock()
        
        with pytest.raises(ValueError, match="Either sample_id or compound_id must be provided"):
            service.create_request(
                db=mock_db,
                requester_id=uuid4(),
                requested_quantity=50.0,
            )

    @patch('amprenta_rag.services.inventory.CompoundRequest')
    def test_create_request_success(self, mock_request):
        """Create request succeeds with valid inputs."""
        mock_db = MagicMock()
        mock_request_obj = MagicMock()
        mock_request_obj.id = uuid4()
        mock_request.return_value = mock_request_obj
        
        result = service.create_request(
            db=mock_db,
            requester_id=uuid4(),
            compound_id=uuid4(),
            requested_quantity=50.0,
            priority="high",
        )
        
        # Verify request was created
        mock_db.add.assert_called_once()
        mock_db.commit.assert_called_once()
        assert result == mock_request_obj

    def test_fulfill_request_validates_quantity(self):
        """P1 TEST: Fulfill with insufficient quantity raises ValueError."""
        mock_db = MagicMock()
        
        # Mock request
        mock_request = MagicMock()
        mock_request.status = "approved"
        mock_request.sample_id = uuid4()
        mock_request.requested_quantity = 100.0
        
        # Mock sample with insufficient quantity
        mock_sample = MagicMock()
        mock_sample.quantity = 50.0  # Less than requested 100
        
        mock_db.query.return_value.filter.return_value.first.side_effect = [
            mock_request,  # Request lookup
            mock_sample,   # Sample lookup
        ]
        
        with pytest.raises(ValueError, match="Insufficient quantity"):
            service.fulfill_request(
                db=mock_db,
                request_id=uuid4(),
                fulfilled_by_id=uuid4(),
            )

    def test_approve_request_validates_status(self):
        """Cannot approve request not in 'requested' status."""
        mock_db = MagicMock()
        
        # Mock already approved request
        mock_request = MagicMock()
        mock_request.status = "approved"  # Not 'requested'
        mock_db.query.return_value.filter.return_value.first.return_value = mock_request
        
        with pytest.raises(ValueError, match="not in 'requested' status"):
            service.approve_request(
                db=mock_db,
                request_id=uuid4(),
                approved_by_id=uuid4(),
            )


class TestAlertQueries:
    """Test alert query functions."""
    
    def test_alert_functions_exist(self):
        """Alert query functions exist and are callable."""
        assert callable(service.get_low_stock_samples)
        assert callable(service.get_expiring_samples) 
        assert callable(service.get_expired_samples)
        
        # Verify function signatures
        import inspect
        
        low_stock_sig = inspect.signature(service.get_low_stock_samples)
        assert 'db' in low_stock_sig.parameters
        assert 'threshold' in low_stock_sig.parameters
        
        expiring_sig = inspect.signature(service.get_expiring_samples)
        assert 'db' in expiring_sig.parameters
        assert 'days' in expiring_sig.parameters


class TestBarcodeLookup:
    """Test unified barcode lookup."""
    
    def test_lookup_barcode_sample_priority(self):
        """Lookup checks samples first, then plates."""
        mock_db = MagicMock()
        
        # Mock sample found
        mock_sample = MagicMock()
        mock_sample.id = uuid4()
        
        with patch.object(service, 'get_compound_sample_by_barcode', return_value=mock_sample):
            with patch.object(service, 'get_plate_by_barcode', return_value=None):
                result = service.lookup_barcode(mock_db, "TEST-BARCODE")
                
                assert result["type"] == "sample"
                assert result["entity"] == mock_sample

    def test_lookup_barcode_plate_fallback(self):
        """Lookup returns plate if no sample found."""
        mock_db = MagicMock()
        
        # Mock plate found, no sample
        mock_plate = MagicMock()
        mock_plate.id = uuid4()
        
        with patch.object(service, 'get_compound_sample_by_barcode', return_value=None):
            with patch.object(service, 'get_plate_by_barcode', return_value=mock_plate):
                result = service.lookup_barcode(mock_db, "TEST-BARCODE")
                
                assert result["type"] == "plate"
                assert result["entity"] == mock_plate

    def test_lookup_barcode_not_found(self):
        """Lookup returns None for unknown barcode."""
        mock_db = MagicMock()
        
        with patch.object(service, 'get_compound_sample_by_barcode', return_value=None):
            with patch.object(service, 'get_plate_by_barcode', return_value=None):
                result = service.lookup_barcode(mock_db, "UNKNOWN-BARCODE")
                
                assert result["type"] is None
                assert result["entity"] is None


class TestServiceFunctionStructure:
    """Test service function structure and behavior."""
    
    def test_create_compound_sample_parameters(self):
        """create_compound_sample has correct parameter structure."""
        # Test that function exists with expected parameters
        import inspect
        sig = inspect.signature(service.create_compound_sample)
        params = list(sig.parameters.keys())
        
        # Verify key parameters exist
        assert 'db' in params
        assert 'compound_id' in params
        assert 'quantity' in params
        assert 'barcode' in params

    def test_fulfill_request_parameters(self):
        """fulfill_request has correct parameter structure."""
        import inspect
        sig = inspect.signature(service.fulfill_request)
        params = list(sig.parameters.keys())
        
        # Verify key parameters exist
        assert 'db' in params
        assert 'request_id' in params
        assert 'fulfilled_by_id' in params

    def test_lookup_barcode_return_structure(self):
        """lookup_barcode returns dict with expected keys."""
        mock_db = MagicMock()
        
        with patch.object(service, 'get_compound_sample_by_barcode', return_value=None):
            with patch.object(service, 'get_plate_by_barcode', return_value=None):
                result = service.lookup_barcode(mock_db, "TEST")
                
                assert isinstance(result, dict)
                assert "type" in result
                assert "entity" in result


class TestServiceValidation:
    """Test service validation logic."""
    
    def test_generate_barcode_timestamp_format(self):
        """Generated barcode includes valid timestamp."""
        barcode = service.generate_barcode()
        parts = barcode.split("-")
        timestamp_part = parts[1]
        
        # Should be valid YYYYMMDD format
        assert len(timestamp_part) == 8
        assert timestamp_part.isdigit()
        
        # Should be reasonable date (within last year to next year)
        year = int(timestamp_part[:4])
        current_year = date.today().year
        assert current_year - 1 <= year <= current_year + 1

    def test_generate_barcode_hex_format(self):
        """Generated barcode includes valid hex suffix."""
        barcode = service.generate_barcode()
        parts = barcode.split("-")
        hex_part = parts[2]
        
        # Should be 6-character uppercase hex
        assert len(hex_part) == 6
        assert hex_part.isupper()
        
        # Should be valid hex
        try:
            int(hex_part, 16)
        except ValueError:
            pytest.fail(f"Invalid hex: {hex_part}")

    def test_service_imports_successfully(self):
        """All service functions import without error."""
        # Verify key functions exist
        assert hasattr(service, 'create_compound_sample')
        assert hasattr(service, 'fulfill_request')
        assert hasattr(service, 'lookup_barcode')
        assert hasattr(service, 'generate_barcode')
        assert hasattr(service, 'get_low_stock_samples')
        assert hasattr(service, 'get_expiring_samples')
        assert hasattr(service, 'create_compound_plate')
        assert hasattr(service, 'create_request')

    def test_service_functions_are_callable(self):
        """All service functions are callable."""
        functions = [
            'create_compound_sample', 'get_compound_sample', 'list_compound_samples',
            'update_compound_sample', 'transfer_compound_sample', 'deplete_compound_sample',
            'create_compound_plate', 'get_plate', 'list_plates', 'get_plate_contents',
            'create_request', 'approve_request', 'fulfill_request', 'reject_request',
            'list_pending_requests', 'get_request', 'list_requests',
            'lookup_barcode', 'get_low_stock_samples', 'get_expiring_samples', 'get_expired_samples'
        ]
        
        for func_name in functions:
            func = getattr(service, func_name)
            assert callable(func), f"{func_name} is not callable"