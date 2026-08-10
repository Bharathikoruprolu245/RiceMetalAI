import uuid
from datetime import datetime

from sqlalchemy import DateTime, String, Text, func
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import Mapped, mapped_column

from backend.database.base import Base


class HeavyMetal(Base):
    __tablename__ = "heavy_metals"

    id: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True),
        primary_key=True,
        default=uuid.uuid4,
    )

    name: Mapped[str] = mapped_column(
        String(100),
        unique=True,
        nullable=False,
    )

    symbol: Mapped[str] = mapped_column(
        String(10),
        unique=True,
        nullable=False,
    )

    ionic_form: Mapped[str | None] = mapped_column(
        String(50),
        nullable=True,
    )

    category: Mapped[str] = mapped_column(
        String(50),
        nullable=False,
    )

    description: Mapped[str | None] = mapped_column(
        Text,
        nullable=True,
    )

    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        server_default=func.now(),
        nullable=False,
    )

    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        server_default=func.now(),
        onupdate=func.now(),
        nullable=False,
    )

    def __repr__(self) -> str:
        return f"<HeavyMetal(name={self.name!r}, symbol={self.symbol!r})>"
