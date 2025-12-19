import jwt
from jupyterhub.auth import Authenticator
from traitlets import Unicode


class TokenAuthenticator(Authenticator):
    """JupyterHub authenticator that validates JWT tokens from Streamlit."""

    jwt_secret = Unicode(
        help="Secret used to validate HS256 JWT tokens"
    ).tag(config=True)

    auto_login = True
    login_service = "Streamlit SSO"

    async def authenticate(self, handler, data):
        """Validate JWT token from query parameter or data dict."""
        # Try query param first (auto_login flow)
        token = handler.get_argument("token", "")

        # Fallback to data dict
        if not token and isinstance(data, dict):
            token = data.get("token", "")

        if not token:
            return None

        try:
            secret = self.jwt_secret or "dev-secret-change-me"
            payload = jwt.decode(token, secret, algorithms=["HS256"])
            username = payload.get("username")
            if username:
                self.log.info(f"Token auth successful for {username}")
                return username
            return None
        except jwt.ExpiredSignatureError:
            self.log.warning("Token expired")
            return None
        except jwt.InvalidTokenError as e:
            self.log.warning(f"Invalid token: {e}")
            return None


