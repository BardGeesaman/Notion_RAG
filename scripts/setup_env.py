#!/usr/bin/env python3
"""Interactive environment setup helper for local development.

This script guides users through setting up their .env file for local
development by copying .env.example and prompting for required values.

Usage:
    python scripts/setup_env.py

Features:
- Interactive prompts for required configuration
- Validates existing .env files
- Provides helpful next steps
- Safe handling of sensitive data
"""

import os
import sys
from pathlib import Path
from typing import Optional


def print_banner():
    """Print the setup banner."""
    print("=" * 60)
    print("🔧 Amprenta RAG - Environment Setup Helper")
    print("=" * 60)
    print()
    print("This helper will guide you through setting up your local")
    print("development environment by configuring your .env file.")
    print()


def print_section(title: str):
    """Print a section header."""
    print(f"\n{'─' * 50}")
    print(f"📋 {title}")
    print("─" * 50)


def get_project_root() -> Path:
    """Get the project root directory."""
    # This script is in scripts/, so project root is parent
    return Path(__file__).parent.parent


def check_files(project_root: Path) -> tuple[Path, Path]:
    """Check for .env.example and .env files."""
    env_example = project_root / ".env.example"
    env_file = project_root / ".env"
    
    if not env_example.exists():
        print(f"❌ Error: .env.example not found at {env_example}")
        print("\nThis file should exist in the project root.")
        print("Please ensure you're running this from the correct directory.")
        sys.exit(1)
    
    print(f"✅ Found .env.example at {env_example}")
    
    return env_example, env_file


def check_existing_env(env_file: Path) -> bool:
    """Check if .env file already exists and get user confirmation."""
    if env_file.exists():
        print(f"\n⚠️  .env file already exists at {env_file}")
        print("\nOptions:")
        print("  y - Overwrite the existing .env file")
        print("  n - Keep existing .env file and exit")
        print("  v - View current .env file first")
        
        while True:
            response = input("\nWhat would you like to do? (y/n/v): ").strip().lower()
            
            if response == 'v':
                print(f"\n📄 Current .env file contents:")
                print("─" * 40)
                try:
                    content = env_file.read_text()
                    # Hide sensitive values for security
                    lines = content.split('\n')
                    for line in lines:
                        if '=' in line and not line.strip().startswith('#'):
                            key, value = line.split('=', 1)
                            if any(secret in key.upper() for secret in ['PASSWORD', 'KEY', 'SECRET', 'TOKEN']):
                                print(f"{key}=***HIDDEN***")
                            else:
                                print(line)
                        else:
                            print(line)
                    print("─" * 40)
                except Exception as e:
                    print(f"Error reading file: {e}")
                continue
            elif response == 'y':
                return True
            elif response == 'n':
                print("\nKeeping existing .env file. Exiting.")
                return False
            else:
                print("Please enter 'y', 'n', or 'v'")
    
    return True


def get_user_input(prompt: str, default: str = "", required: bool = False, 
                  secret: bool = False, validate_fn=None) -> str:
    """Get user input with validation."""
    while True:
        if default:
            full_prompt = f"{prompt} [{default}]: "
        else:
            full_prompt = f"{prompt}: "
        
        if secret:
            import getpass
            value = getpass.getpass(full_prompt)
        else:
            value = input(full_prompt).strip()
        
        if not value and default:
            value = default
        
        if required and not value:
            print("❌ This field is required. Please enter a value.")
            continue
        
        if validate_fn and value:
            if not validate_fn(value):
                continue
        
        return value


def validate_email(email: str) -> bool:
    """Basic email validation."""
    if '@' not in email or '.' not in email:
        print("❌ Please enter a valid email address.")
        return False
    return True


def validate_openai_key(key: str) -> bool:
    """Basic OpenAI API key validation."""
    if not key.startswith(('sk-', 'sk-proj-')):
        print("❌ OpenAI API keys should start with 'sk-' or 'sk-proj-'")
        return False
    if len(key) < 20:
        print("❌ OpenAI API key seems too short.")
        return False
    return True


def setup_database_config() -> dict[str, str]:
    """Setup database configuration."""
    print_section("Database Configuration")
    print("📊 PostgreSQL database settings (required for core functionality)")
    print()
    
    config = {}
    
    # Database connection details
    config['POSTGRES_USER'] = get_user_input(
        "Database username", 
        default="postgres"
    )
    
    config['POSTGRES_PASSWORD'] = get_user_input(
        "Database password", 
        required=True,
        secret=True
    )
    
    config['POSTGRES_HOST'] = get_user_input(
        "Database host", 
        default="localhost"
    )
    
    config['POSTGRES_PORT'] = get_user_input(
        "Database port", 
        default="5432"
    )
    
    config['POSTGRES_DB'] = get_user_input(
        "Database name", 
        default="amprenta"
    )
    
    return config


def setup_api_keys() -> dict[str, str]:
    """Setup API keys for external services."""
    print_section("API Keys (Optional)")
    print("🔑 Configure API keys for external service integrations")
    print("   (You can skip these and add them later)")
    print()
    
    config = {}
    
    # OpenAI API Key
    print("🤖 OpenAI API Key (for LLM features like RAG, chat, summarization)")
    openai_key = get_user_input(
        "OpenAI API key (press Enter to skip)",
        validate_fn=lambda k: validate_openai_key(k) if k else True
    )
    if openai_key:
        config['OPENAI_API_KEY'] = openai_key
    
    # NCBI Email
    print("\n📧 NCBI Email (required by NCBI for API access to GEO, PubMed)")
    ncbi_email = get_user_input(
        "NCBI email address (press Enter to skip)",
        validate_fn=lambda e: validate_email(e) if e else True
    )
    if ncbi_email:
        config['NCBI_EMAIL'] = ncbi_email
    
    # Zotero API Key
    print("\n📚 Zotero API Key (for bibliography and reference management)")
    print("   Get your key from: https://www.zotero.org/settings/keys")
    zotero_key = get_user_input(
        "Zotero API key (press Enter to skip)"
    )
    if zotero_key:
        config['ZOTERO_API_KEY'] = zotero_key
        
        zotero_library = get_user_input(
            "Zotero library ID (press Enter to skip)"
        )
        if zotero_library:
            config['ZOTERO_LIBRARY_ID'] = zotero_library
    
    # GEO API Key
    print("\n🧬 GEO API Key (for higher rate limits on NCBI GEO queries)")
    geo_key = get_user_input(
        "GEO API key (press Enter to skip)"
    )
    if geo_key:
        config['GEO_API_KEY'] = geo_key
    
    return config


def setup_development_options() -> dict[str, str]:
    """Setup development-specific options."""
    print_section("Development Options")
    print("⚙️  Configure development-specific settings")
    print()
    
    config = {}
    
    # Environment
    config['ENVIRONMENT'] = get_user_input(
        "Environment designation",
        default="dev"
    )
    
    # Debug mode
    debug_mode = get_user_input(
        "Enable debug mode? (y/n)",
        default="n"
    ).lower() == 'y'
    config['DEBUG'] = "true" if debug_mode else "false"
    
    # Log level
    if debug_mode:
        config['LOG_LEVEL'] = "DEBUG"
    else:
        log_level = get_user_input(
            "Log level (DEBUG/INFO/WARNING/ERROR)",
            default="INFO"
        ).upper()
        config['LOG_LEVEL'] = log_level
    
    # Disable auth for local dev
    disable_auth = get_user_input(
        "Disable authentication for local development? (y/n)",
        default="y"
    ).lower() == 'y'
    config['DISABLE_AUTH'] = "true" if disable_auth else "false"
    
    return config


def update_env_content(content: str, config: dict[str, str]) -> str:
    """Update .env.example content with user-provided values."""
    lines = content.split('\n')
    updated_lines = []
    
    for line in lines:
        if '=' in line and not line.strip().startswith('#'):
            key = line.split('=', 1)[0]
            if key in config:
                updated_lines.append(f"{key}={config[key]}")
            else:
                updated_lines.append(line)
        else:
            updated_lines.append(line)
    
    return '\n'.join(updated_lines)


def write_env_file(env_file: Path, content: str):
    """Write the .env file safely."""
    try:
        env_file.write_text(content)
        print(f"\n✅ Successfully created {env_file}")
    except Exception as e:
        print(f"\n❌ Error writing .env file: {e}")
        sys.exit(1)


def print_next_steps(project_root: Path):
    """Print helpful next steps for the user."""
    print_section("Next Steps")
    print("🚀 Your .env file is ready! Here's what to do next:")
    print()
    
    steps = [
        "Review and edit .env file if needed:",
        f"   nano {project_root}/.env",
        "",
        "Activate your conda environment:",
        "   conda activate myenv",
        "",
        "Start PostgreSQL database:",
        "   docker-compose up -d postgres",
        "   # OR install PostgreSQL locally",
        "",
        "Run database migrations:",
        "   alembic upgrade head",
        "",
        "Start the API server:",
        "   python -m uvicorn amprenta_rag.api.main:app --reload",
        "",
        "Start the dashboard (in another terminal):",
        "   streamlit run scripts/dashboard/main.py",
        "",
        "Test your setup:",
        "   curl http://localhost:8000/health",
        "   # Should return: {\"status\": \"healthy\"}"
    ]
    
    for i, step in enumerate(steps, 1):
        if step:
            if step.startswith("   "):
                print(f"     {step[3:]}")
            else:
                print(f"{i:2d}. {step}")
        else:
            print()


def print_troubleshooting():
    """Print troubleshooting information."""
    print_section("Troubleshooting")
    print("🔧 If you encounter issues:")
    print()
    print("• Database connection errors:")
    print("    Check PostgreSQL is running: docker ps")
    print("    Verify credentials in .env file")
    print()
    print("• API key errors:")
    print("    Test OpenAI key: python -c \"from amprenta_rag.utils.secrets import get_api_key; print(get_api_key('openai'))\"")
    print("    Check key format and permissions")
    print()
    print("• Import errors:")
    print("    Ensure conda environment is activated")
    print("    Install dependencies: pip install -r requirements.txt")
    print()
    print("• For more help, see docs/SECRETS_MANAGEMENT.md")


def main():
    """Main setup function."""
    try:
        print_banner()
        
        # Check file locations
        project_root = get_project_root()
        env_example, env_file = check_files(project_root)
        
        # Check existing .env file
        if not check_existing_env(env_file):
            return
        
        # Read .env.example template
        try:
            content = env_example.read_text()
        except Exception as e:
            print(f"❌ Error reading .env.example: {e}")
            sys.exit(1)
        
        # Collect configuration
        all_config = {}
        
        # Database configuration (required)
        db_config = setup_database_config()
        all_config.update(db_config)
        
        # API keys (optional)
        api_config = setup_api_keys()
        all_config.update(api_config)
        
        # Development options (optional)
        dev_config = setup_development_options()
        all_config.update(dev_config)
        
        # Update content with user values
        updated_content = update_env_content(content, all_config)
        
        # Write .env file
        write_env_file(env_file, updated_content)
        
        # Print next steps
        print_next_steps(project_root)
        print_troubleshooting()
        
        print("\n" + "=" * 60)
        print("🎉 Environment setup complete!")
        print("=" * 60)
        
    except KeyboardInterrupt:
        print("\n\n⚠️  Setup cancelled by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
