from logging.config import fileConfig

from alembic import context
from sqlalchemy import pool

from backend.database.base import Base
from backend.database.session import engine

# Import models so Alembic can detect them.
# More models will be added here as the project grows.
  # noqa: F401
from backend.models import (
    metal,
    gene_family,
    gene,
    gene_metal,
)

config = context.config

if config.config_file_name is not None:
    fileConfig(config.config_file_name)

# Alembic uses this metadata for autogenerate.
target_metadata = Base.metadata


def run_migrations_offline() -> None:
    """Run migrations without creating a database connection."""

    url = str(engine.url)

    context.configure(
        url=url,
        target_metadata=target_metadata,
        literal_binds=True,
        dialect_opts={"paramstyle": "named"},
    )

    with context.begin_transaction():
        context.run_migrations()


def run_migrations_online() -> None:
    """Run migrations using the configured SQLAlchemy engine."""

    with engine.connect() as connection:
        context.configure(
            connection=connection,
            target_metadata=target_metadata,
        )

        with context.begin_transaction():
            context.run_migrations()


if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()
