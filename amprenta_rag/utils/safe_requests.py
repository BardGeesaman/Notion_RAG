"""Safe HTTP request utilities with SSRF protection."""

import ipaddress
import socket
from typing import Optional
from urllib.parse import urlparse

import requests

# Blocked schemes
BLOCKED_SCHEMES = {"file", "ftp", "gopher", "data", "javascript"}

# Private IP ranges
PRIVATE_IP_RANGES = [
    ipaddress.ip_network("10.0.0.0/8"),
    ipaddress.ip_network("172.16.0.0/12"),
    ipaddress.ip_network("192.168.0.0/16"),
    ipaddress.ip_network("127.0.0.0/8"),
    ipaddress.ip_network("169.254.0.0/16"),
    ipaddress.ip_network("0.0.0.0/8"),
    ipaddress.ip_network("::1/128"),
    ipaddress.ip_network("fe80::/10"),
    ipaddress.ip_network("fc00::/7"),
]

# Allowed external domains (whitelist for known APIs)
ALLOWED_DOMAINS = {
    "api.ncbi.nlm.nih.gov",
    "eutils.ncbi.nlm.nih.gov",
    "pubchem.ncbi.nlm.nih.gov",
    "www.ebi.ac.uk",
    "rest.ensembl.org",
    "rest.kegg.jp",
    "api.openalex.org",
    "api.semanticscholar.org",
    "alphafold.ebi.ac.uk",
    "files.rcsb.org",
    "string-db.org",
    "www.uniprot.org",
}


def is_safe_url(url: str, allow_internal: bool = False) -> tuple[bool, str]:
    """
    Validate URL is safe for external requests.
    
    Args:
        url: URL to validate
        allow_internal: If True, skip private IP checks (for internal services)
    
    Returns:
        Tuple of (is_safe, error_message)
    """
    try:
        parsed = urlparse(url)
    except Exception:
        return False, "Invalid URL format"
    
    # Block dangerous schemes
    if parsed.scheme.lower() in BLOCKED_SCHEMES:
        return False, f"Blocked scheme: {parsed.scheme}"
    
    hostname = parsed.hostname
    if not hostname:
        return False, "Invalid URL: no hostname"
    
    # Block localhost
    if hostname.lower() in ("localhost", "0.0.0.0"):
        return False, "Blocked: localhost"
    
    if not allow_internal:
        # Check private IPs
        try:
            ip = ipaddress.ip_address(hostname)
            for network in PRIVATE_IP_RANGES:
                if ip in network:
                    return False, "Blocked: private IP range"
        except ValueError:
            # Hostname, resolve and check
            try:
                ip = ipaddress.ip_address(socket.gethostbyname(hostname))
                for network in PRIVATE_IP_RANGES:
                    if ip in network:
                        return False, "Blocked: private IP range"
            except (socket.gaierror, ValueError):
                pass
    
    return True, ""


def safe_get(url: str, allow_internal: bool = False, **kwargs) -> requests.Response:
    """
    Safe GET request with SSRF protection.
    
    Raises:
        ValueError: If URL is blocked
    """
    is_safe, error = is_safe_url(url, allow_internal)
    if not is_safe:
        raise ValueError(f"URL blocked: {error}")
    
    return requests.get(url, **kwargs)


def safe_post(url: str, allow_internal: bool = False, **kwargs) -> requests.Response:
    """Safe POST request with SSRF protection."""
    is_safe, error = is_safe_url(url, allow_internal)
    if not is_safe:
        raise ValueError(f"URL blocked: {error}")
    
    return requests.post(url, **kwargs)
